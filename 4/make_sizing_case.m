function [sizing, P_sized, S_worst] = make_sizing_case(scaled, params, gainsPF, opts)
%% DETERMÄ°NÄ°STÄ°K BOYUTLANDIRMA
% GA sonrasÄ±nda elde edilen kazanÃ§lar ve PF parametreleri sabitlenerek
% damper orifisinin boyutlandÄ±rÄ±lmasÄ± yapÄ±lÄ±r. Hesaplamalar, GA sÄ±rasÄ±nda
% kullanÄ±lan ODE modeline uyumlu basit kurallara dayanÄ±r.
%
% Girdiler
%   scaled  - GA'da kullanÄ±lan Ã¶lÃ§eklenmiÅŸ yer hareketleri
%   params  - yapÄ± + damper temel parametre yapÄ±sÄ±
%   gainsPF - {g_lo,g_mid,g_hi, PF_tau, PF_gain} alanlarÄ±nÄ± iÃ§eren yapÄ±
%   opts    - isteÄŸe baÄŸlÄ± ayarlar:
%               dp_allow_frac (varsayÄ±lan 0.60)
%               alpha_lam     (varsayÄ±lan 0.15)
%               n_orf_set     (varsayÄ±lan [5 6])
%               dmm_step      (varsayÄ±lan 0.05 mm)
%               Cd_init       (varsayÄ±lan params.orf.CdInf ya da 0.70)
%               mu_nom        (varsayÄ±lan params.mu_ref)
%
% Ã‡Ä±ktÄ±lar
%   sizing  - tÃ¼retilen bÃ¼yÃ¼klÃ¼kler (Q95_worst, Ao, d_o, n_orf, Cd, Lori, ...)
%   P_sized - {n_orf, orf.d_o, A_o, Qcap_big} alanlarÄ± gÃ¼ncellenmiÅŸ parametre
%             yapÄ±sÄ±, kazanÃ§lar/PF sabitlenmiÅŸ
%   S_worst - run_batch_windowed(scaled, P_sized, O) sonucu rapor yapÄ±sÄ±
%
% -------------------------------------------------------------------------

    % 0 argÃ¼man/az argÃ¼man kolaylÄ±ÄŸÄ±: Ã§alÄ±ÅŸma alanÄ± ve son GA CSV'lerinden verileri Ã§ekmeyi dene
    if nargin < 3 || isempty(scaled) || isempty(params) || isempty(gainsPF)
        % scaled / params deÄŸiÅŸkenlerini temel Ã§alÄ±ÅŸma alanÄ±ndan al
        try
            if nargin < 1 || isempty(scaled), scaled = evalin('base','scaled'); end
        catch
        end
        try
            if nargin < 2 || isempty(params), params = evalin('base','params'); end
        catch
        end
        % En iyi kazanÃ§lar/PF iÃ§in son ga_knee.csv veya ga_front.csv dosyasÄ±nÄ± kullan
        if nargin < 3 || isempty(gainsPF)
            try
                dd = dir(fullfile('out','ga_*'));
                assert(~isempty(dd), 'make_sizing_case:auto: no GA output under out/ga_*');
                [~,ix] = max([dd.datenum]);
                latestDir = fullfile(dd(ix).folder, dd(ix).name);
                kneeCsv = fullfile(latestDir,'ga_knee.csv');
                if exist(kneeCsv,'file')
                    Tbest = readtable(kneeCsv);
                    row = 2; if height(Tbest) < 2, row = 1; end  % baseline at 1, knee at 2
                else
                    frontCsv = fullfile(latestDir,'ga_front.csv');
                    assert(exist(frontCsv,'file')==2, 'make_sizing_case:auto: ga_front.csv not found');
                    Tbest = readtable(frontCsv);
                    % knee selection on-the-fly
                    f1v = Tbest.f1; f2v = Tbest.f2;
                    f1n = (f1v - min(f1v)) ./ max(eps, (max(f1v)-min(f1v)));
                    f2n = (f2v - min(f2v)) ./ max(eps, (max(f2v)-min(f2v)));
                    [~,row] = min(hypot(f1n, f2n));
                end
                req = {'g_lo','g_mid','g_hi','PF_tau','PF_gain'};
                assert(all(ismember(req, Tbest.Properties.VariableNames)), ...
                    'make_sizing_case:auto: required columns missing in GA CSV');
                gainsPF = struct('g_lo',Tbest.g_lo(row), 'g_mid',Tbest.g_mid(row), ...
                                 'g_hi',Tbest.g_hi(row), 'PF_tau',Tbest.PF_tau(row), ...
                                 'PF_gain',Tbest.PF_gain(row));
            catch ME
                error('make_sizing_case:auto_init','Failed to auto-infer gains/PF from GA outputs: %s', ME.message);
            end
        end
    end

    if nargin < 4 || isempty(opts), opts = struct(); end

    % --- varsayÄ±lanlÄ± pratik getter ---
    function v = getopt(s, f, def)
        if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
            v = s.(f);
        else
            v = def;
        end
    end

    %% VarsayÄ±lan AyarlarÄ±n Ã‡ekilmesi
    dp_allow_frac = getopt(opts,'dp_allow_frac',0.60);   % Î”p_kv,target = 0.6*Î”p_cap
    alpha_lam     = getopt(opts,'alpha_lam',0.15);       % laminer/tÃ¼rbÃ¼lans oran hedefi
    n_orf_set     = getopt(opts,'n_orf_set',[5 6]);      % izinli delik sayÄ±larÄ±
    dmm_step      = getopt(opts,'dmm_step',0.05);        % d_o kuantizasyon adÄ±mÄ± (mm)
    Cd_init       = getopt(opts,'Cd_init',getopt(getopt(params,'orf',struct()),'CdInf',0.70));
    mu_nom        = getopt(opts,'mu_nom',getopt(params,'mu_ref',0.9));  % PaÂ·s

    %% AdÄ±m 1: PF ve Geometri KazanÃ§larÄ±nÄ± Sabitleme
    P = params;
    % kazanÃ§lar -> toggle_gain (alt 3, orta, Ã¼st 2 kat)
    nStories = size(P.M,1)-1;
    tg = ones(nStories,1) * gainsPF.g_mid;
    loN = min(3,nStories); if loN>0, tg(1:loN) = gainsPF.g_lo; end
    hiN = min(2,nStories); if hiN>0, tg(end-hiN+1:end) = gainsPF.g_hi; end
    P.toggle_gain = tg;
    % PF parametreleri
    if isfield(P,'cfg') && isfield(P.cfg,'PF')
        P.cfg.PF.tau  = gainsPF.PF_tau;
        P.cfg.PF.gain = gainsPF.PF_gain;
    end

    %% AdÄ±m 2: Sabit KazanÃ§larla Sistem DeÄŸerlendirmesi
    O = struct('do_export',false,'quiet',true,'thermal_reset','each','order','natural', ...
               'use_orifice',true,'use_thermal',true, ...
               'mu_factors',[0.75 1.00 1.25], 'mu_weights',[0.2 0.6 0.2], 'thr', []);
    S_worst = run_batch_windowed(scaled, P, O);

    %% AdÄ±m 3: Q95 ve BasÄ±nÃ§ SÄ±nÄ±rlarÄ±nÄ±n Belirlenmesi
    [Q95_worst, dp_kv_target, dp_cap] = find_Q95_worst(S_worst, O, dp_allow_frac);

    %% AdÄ±m 4: Ao/d_o AramasÄ± ve Laminer Kol
    rho = P.rho;
    best = select_orifice_size(Q95_worst, dp_kv_target, Cd_init, n_orf_set, dmm_step, alpha_lam, mu_nom, rho);

    %% AdÄ±m 5: Reynolds SayÄ±sÄ± ve Cd GÃ¼ncellemesi
    best = update_Re_Cd(best, Q95_worst, dp_kv_target, mu_nom, P, dmm_step, rho);

    %% AdÄ±m 6: Debi Ãœst SÄ±nÄ±rÄ±nÄ±n Belirlenmesi
    Qcap_big = 1.2 * abs(Q95_worst);

    %% AdÄ±m 7: Parametre YapÄ±sÄ±nÄ±n GÃ¼ncellenmesi
    P_sized = P;
    P_sized.n_orf    = best.n_orf;
    P_sized.orf.d_o  = best.d_o;
    P_sized.A_o      = best.A_o;
    P_sized.Qcap_big = Qcap_big;

    %% AdÄ±m 8: BoyutlandÄ±rma Paketinin Derlenmesi
    sizing = struct();
    sizing.Q95_worst   = Q95_worst;
    sizing.dp_cap      = dp_cap;
    sizing.dp_kv_target= dp_kv_target;
    sizing.n_orf       = best.n_orf;
    sizing.d_o         = best.d_o;
    sizing.A_o         = best.A_o;
    sizing.Cd          = best.Cd;
    sizing.L_orif      = best.Lori;
    sizing.Qcap_big    = Qcap_big;
    sizing.alpha_lam   = alpha_lam;
    sizing.notes       = 'Deterministic sizing pass; PF/gains fixed.';

    %% AdÄ±m 9: Parametre FarklarÄ±nÄ±n RaporlanmasÄ±
    try
        [updates, T] = sizing_param_diff(params, P_sized, gainsPF, sizing); %#ok<NASGU>
        % Otomatik yamaya hazÄ±r satÄ±rlar
        fprintf('parametre.m iÃ§in Ã¶nerilen satÄ±rlar:\n');
        for i=1:numel(updates.lines)
            fprintf('  %s\n', updates.lines{i});
        end\n        try\n            report_toggle = getopt(opts,'report_toggle',false);\n        catch\n            report_toggle = false;\n        end\n        if report_toggle\n            fprintf(''  %% toggle_gain vekter:\n  toggle_gain = %s;  %% (n-1)x1\n'', mat2str(updates.toggle_gain_vec,3));\n        end\n
    catch ME
        warning('sizing_param_diff failed: %s', ME.message);
    end
end

%% Yerel fonksiyon: Q95_worst ve basÄ±nÃ§ sÄ±nÄ±rÄ±nÄ± bulur
function [Q95_worst, dp_kv_target, dp_cap] = find_Q95_worst(S_worst, O, dp_allow_frac)
    % Ã–zet tablo Ã¼zerinden en kÃ¶tÃ¼ %95 debiyi Ã§ek
    Q95_worst = NaN;
    try
        vars = S_worst.table.Properties.VariableNames;
        if ismember('Q_q95_worst', vars)
            v = S_worst.table.Q_q95_worst;
            Q95_worst = max(v(:));
        elseif ismember('Q_q95_w', vars)
            v = S_worst.table.Q_q95_w;
            Q95_worst = max(v(:));
        end
    catch
    end
    % Model hard-limit veya thr Ã¼zerinden kap
    dp_cap_default = 1.0e9; % Pa
    dp_cap = dp_cap_default;
    if isfield(O,'thr') && ~isempty(O.thr) && isfield(O.thr,'dP95_max') && ~isempty(O.thr.dP95_max)
        dp_cap = O.thr.dP95_max;
    end
    dp_kv_target = dp_allow_frac * dp_cap;
end

%% Yerel fonksiyon: Ao/d_o seÃ§imi ve laminer kol hesabÄ±
function best = select_orifice_size(Q95_worst, dp_kv_target, Cd_init, n_orf_set, dmm_step, alpha_lam, mu_nom, rho)
    % TÃ¼rbÃ¼lans baskÄ±n koÅŸuldan gerekli alan
    Ao_req = abs(Q95_worst) / max(Cd_init,eps) * sqrt(rho / max(2*dp_kv_target,eps));
    dgrid = @(dmm) (dmm_step * round(dmm./max(dmm_step,eps))); % mm Ä±zgara
    best = struct('n_orf',NaN,'d_o',NaN,'A_o',NaN,'Cd',Cd_init,'Lori',NaN);
    best_err = inf;
    for n_orf = n_orf_set(:).' % izinli delik sayÄ±larÄ± Ã¼zerinde ara
        d_o = sqrt(4*Ao_req/(pi*n_orf));       % m
        dmm = dgrid(1000*d_o); d_o_q = dmm/1000;  % kuantize m
        A_o = n_orf * (pi*d_o_q^2/4);
        err = abs(A_o - Ao_req)/max(Ao_req,eps);
        if err < best_err
            best_err = err;
            best.n_orf = n_orf; best.d_o = d_o_q; best.A_o = A_o;
        end
    end
    % Laminer kol uzunluÄŸu
    L_orif = alpha_lam * dp_kv_target * (pi * best.d_o^4) * best.n_orf / max(128*mu_nom*abs(Q95_worst),eps);
    best.Lori = L_orif;
end

%% Yerel fonksiyon: Reynolds sayÄ±sÄ± ve Cd gÃ¼ncellemesi
function best = update_Re_Cd(best, Q95_worst, dp_kv_target, mu_nom, P, dmm_step, rho)
    dgrid = @(dmm) (dmm_step * round(dmm./max(dmm_step,eps))); % mm Ä±zgara
    Cd = best.Cd; Re = NaN; %#ok<NASGU>
    for it = 1:2 % tekrarlÄ± gÃ¼ncelleme
        Q_hole = Q95_worst / max(best.n_orf,1);
        Re = 4*rho*abs(Q_hole) / max(pi*mu_nom*best.d_o,eps);
        Cd = P.orf.CdInf - (P.orf.CdInf - P.orf.Cd0) / (1 + (Re/max(P.orf.Rec,eps))^P.orf.p_exp);
        Ao_req = abs(Q95_worst) / max(Cd,eps) * sqrt(rho / max(2*dp_kv_target,eps));
        d_o = sqrt(4*Ao_req/(pi*best.n_orf));
        dmm = dgrid(1000*d_o); best.d_o = dmm/1000; best.A_o = best.n_orf*(pi*best.d_o^2/4);
    end
    best.Cd = Cd;
end

%% ============================================================
% ALT FONKSÄ°YON: sizing_param_diff
% Parametre dosyasÄ±ndaki deÄŸiÅŸiklikleri Ã¶zetler ve raporlar.
%% ============================================================
function [updates, T] = sizing_param_diff(P_old, P_sized, gainsPF, sizing)
% Raporda sadece deÄŸiÅŸen parametreler (eski -> yeni) gÃ¶sterilir.

    % --- gÃ¼venli getter -------------------------------------------------
    function v = getf(s,f,def)
        if nargin<3, def = []; end
        if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
            v = s.(f);
        else
            v = def;
        end
    end

    % --- Eski deÄŸerleri topla -------------------------------------------
    old_d_o   = getf(P_old,'d_o',  getf(getf(P_old,'orf',struct()),'d_o',NaN));
    old_Lori  = getf(P_old,'Lori', getf(getf(P_old,'orf',struct()),'L_orif',NaN));
    old_n_orf = getf(P_old,'n_orf',NaN);
    % Eksik GA anlÄ±k gÃ¶rÃ¼ntÃ¼lerinde Lori ve n_orf tahmini
    if ~isfinite(old_n_orf)
        old_Ao_tmp = getf(P_old,'A_o',NaN);
        if isfinite(old_Ao_tmp) && isfinite(old_d_o)
            old_n_orf = old_Ao_tmp / (pi*(old_d_o^2)/4);
        end
    end
    if ~isfinite(old_Lori)
        c_lam0 = getf(P_old,'c_lam0',NaN);
        mu_ref = getf(P_old,'mu_ref',NaN);
        Ap     = getf(P_old,'Ap',NaN);
        if all(isfinite([c_lam0 mu_ref Ap old_d_o]))
            old_Lori = c_lam0*(old_d_o^4)/(12*mu_ref*(Ap^2));
        end
    end
    if isfinite(old_n_orf) && isfinite(old_d_o)
        old_Ao = old_n_orf*pi*(old_d_o^2)/4;
    else
        old_Ao = getf(P_old,'A_o',NaN);
    end
    old_Qcap  = getf(P_old,'Qcap_big',NaN);
    old_tau   = getf(getf(getf(P_old,'cfg',struct()),'PF',struct()),'tau',NaN);
    old_gain  = getf(getf(getf(P_old,'cfg',struct()),'PF',struct()),'gain',NaN);
    old_tg    = getf(P_old,'toggle_gain',[]);

    % --- Yeni/Ã¶nerilen deÄŸerleri topla ---------------------------------
    new_d_o   = getf(getf(P_sized,'orf',struct()),'d_o', getf(P_sized,'d_o',NaN));
    new_Lori  = getf(getf(P_sized,'orf',struct()),'L_orif', getf(P_sized,'Lori', getf(sizing,'L_orif',NaN)));
    new_n_orf = getf(P_sized,'n_orf',NaN);
    if isfinite(new_n_orf) && isfinite(new_d_o)
        new_Ao = new_n_orf*pi*(new_d_o^2)/4;
    else
        new_Ao = getf(P_sized,'A_o',NaN);
    end
    new_Qcap  = getf(P_sized,'Qcap_big',NaN);
    new_tau   = getf(gainsPF,'PF_tau',NaN);
    new_gain  = getf(gainsPF,'PF_gain',NaN);

    % --- toggle_gain vektÃ¶rÃ¼ (alt=3, Ã¼st=2 kat) -------------------------
    try
        nStories = size(P_old.M,1)-1;
    catch
        nStories = numel(old_tg); if isempty(nStories) || nStories==0, nStories = 10; end
    end
    old_tg = old_tg(:);
    if isempty(old_tg)
        old_tg = ones(nStories,1)*getf(gainsPF,'g_mid',1);
    elseif isscalar(old_tg)
        old_tg = repmat(old_tg,nStories,1);
    end
    tg = ones(nStories,1)*getf(gainsPF,'g_mid',1);
    loN = min(3,nStories); if loN>0, tg(1:loN) = getf(gainsPF,'g_lo',1); end
    hiN = min(2,nStories); if hiN>0, tg(end-hiN+1:end) = getf(gainsPF,'g_hi',1); end
    new_tg = tg;

    % --- Skaler eÅŸitlik denetimi ---------------------------------------
    function tf = eqnum(a,b)
        if ~(isnumeric(a) && isnumeric(b)) || isempty(a) || isempty(b)
            tf = false; return;
        end
        a1 = a(1); b1 = b(1); % skalerleÅŸtir
        if ~(isfinite(a1) && isfinite(b1))
            tf = false; return;
        end
        tol = max(1e-12, 1e-6*max(1, max(abs([a1 b1]))));
        tf  = abs(a1-b1) <= tol;
    end

    changed = @(oldv,newv) ~( (isnumeric(oldv)&&isnumeric(newv)&&eqnum(oldv,newv)) || isequal(oldv,newv) );

    rows = {
     'd_o (m)',           old_d_o,   new_d_o,   'd_o = %.6g;';
     'Lori (m)',          old_Lori,  new_Lori,  'Lori = %.6g;';
     'n_orf (-)',         old_n_orf, new_n_orf, 'n_orf = %d;';
     'A_o (m^2)',         old_Ao,    new_Ao,    'A_o = %.6g;';
     'Qcap_big (m^3/s)',  old_Qcap,  new_Qcap,  'Qcap_big = %.6g;';
     'cfg.PF.tau (s)',    old_tau,   new_tau,   'cfg.PF.tau = %.6g;';
     'cfg.PF.gain (-)',   old_gain,  new_gain,  'cfg.PF.gain = %.6g'
    };

    Name = {}; Old = {}; New = {}; Template = {};
    for i=1:size(rows,1)
        if changed(rows{i,2}, rows{i,3})
            Name{end+1,1} = rows{i,1}; %#ok<AGROW>
            if contains(rows{i,1},"n_orf")
                Old{end+1,1} = sprintf('%d', round(rows{i,2})); %#ok<AGROW>
                New{end+1,1} = sprintf('%d', round(rows{i,3})); %#ok<AGROW>
            else
                Old{end+1,1} = sprintf('%.6g', rows{i,2}); %#ok<AGROW>
                New{end+1,1} = sprintf('%.6g', rows{i,3}); %#ok<AGROW>
            end
            Template{end+1,1} = rows{i,4}; %#ok<AGROW>
        end
    end
    T = table(Name, Old, New, Template);

    % --- toggle_gain farkÄ± ----------------------------------------------
    tg_changed = numel(old_tg)~=numel(new_tg) || any(abs(old_tg(:)-new_tg(:))>1e-6);
    if tg_changed
        T = [T; table({"toggle_gain ((n-1)x1)"}, {mat2str(old_tg,3)}, {mat2str(new_tg,3)}, {''}, 'VariableNames', T.Properties.VariableNames)];
        fprintf('\n[update] toggle_gain:\n  old: %s\n  new: %s\n', mat2str(old_tg,3), mat2str(new_tg,3));
    end

    % --- konsol Ã§Ä±ktÄ±sÄ± -------------------------------------------------
    fprintf('\n=== PARAMETRE GÃœNCELLEME Ã–NERÄ°LERÄ° (sadece deÄŸiÅŸenler) ===\n');
    for i=1:height(T)
        nm = string(T.Name{i});
        if T.Template{i}~=""
            if contains(nm,"n_orf")
                fprintf(' - %-18s : %s  -->  %s   (parametre.m: %s)\n', nm, T.Old{i}, T.New{i}, sprintf(T.Template{i}, round(str2double(T.New{i}))));
            else
                fprintf(' - %-18s : %s  -->  %s   (parametre.m: %s)\n', nm, T.Old{i}, T.New{i}, sprintf(T.Template{i}, str2double(T.New{i})));
            end
        else
            fprintf(' - %-18s : %s  -->  %s   (vektÃ¶r)\n', nm, T.Old{i}, T.New{i});
        end
    end
    fprintf('===========================================================\n\n');

    % --- CSV yazÄ±mÄ± -----------------------------------------------------
    outdir = fullfile('out'); if ~exist(outdir,'dir'), mkdir(outdir); end
    safe_write(T, fullfile(outdir,'sizing_updates.csv'), @writetable);

    % --- parametre.m'e yapÄ±ÅŸtÄ±rmalÄ±k satÄ±rlar --------------------------
    updates = struct(); updates.lines = {};
    for i=1:height(T)
        if T.Template{i}~=""
            if contains(T.Name{i},"n_orf")
                updates.lines{end+1} = sprintf(T.Template{i}, round(str2double(T.New{i}))); %#ok<AGROW>
            else
                updates.lines{end+1} = sprintf(T.Template{i}, str2double(T.New{i})); %#ok<AGROW>
            end
        end
    end
    updates.toggle_gain_vec = new_tg;
end

function safe_write(obj, filepath, writeFcn)
    try
        writeFcn(obj, filepath);
    catch ME
        warning('Yazma hatasÄ± (%s): %s', filepath, ME.message);
    end
end


