function [sizing, P_sized, S_worst] = make_sizing_case(scaled, params, gainsPF, opts)
%MAKE_SIZING_CASE GA sonrası sabitlenen kazancı ve PF ayarlarını
% Not: Bu sürüm ana dal ile çakışmalar çözülerek güncellenmiştir.
% kullanarak deterministik boyutlandırma yapar ve parametre farklarını
% raporlar.
%
% [SIZING, P_SIZED, S_WORST] = MAKE_SIZING_CASE(SCALED, PARAMS, GAINSPF, OPTS)
%   SCALED yer hareketi seti, PARAMS temel parametreleri ve GAINSPF
%   kazanc/PF değerlerini kullanır. OPS opsiyonel ayarları içerir.
%
% GAINSPF alanları: g_lo, g_mid, g_hi, PF_tau, PF_gain
% OPTS alanları (opsiyonel):
%   dp_allow_frac - Δp_kv,target = dp_allow_frac * Δp_cap (vars. 0.60)
%   alpha_lam     - laminer/türbülans oranı (vars. 0.15)
%   n_orf_set     - delik sayısı adayları (vars. [5 6])
%   dmm_step      - d_o için mm kılavuzu (vars. 0.05)
%   Cd_init       - başlangıç Cd (vars. params.orf.CdInf)
%   mu_nom        - nominal viskozite [Pa·s] (vars. params.mu_ref)

%% ------------------------------------------------------------------------
%% 1) GİRİŞ KONTROLLERİ VE OTOMATİK ÇEKME
if nargin < 3 || isempty(scaled) || isempty(params) || isempty(gainsPF)
    % scaled / params çalışma alanından alınır
    try, if nargin < 1 || isempty(scaled),  scaled  = evalin('base','scaled');  end; end
    try, if nargin < 2 || isempty(params),  params  = evalin('base','params');  end; end
    % GA sonu çıktılarından en iyi kazanc/PF değerleri
    if nargin < 3 || isempty(gainsPF)
        dd = dir(fullfile('out','ga_*'));
        assert(~isempty(dd), 'GA çıktısı bulunamadı.');
        [~,ix] = max([dd.datenum]);
        latestDir = fullfile(dd(ix).folder, dd(ix).name);
        kneeCsv = fullfile(latestDir,'ga_knee.csv');
        if exist(kneeCsv,'file')
            Tbest = readtable(kneeCsv);
            row = 2; if height(Tbest) < 2, row = 1; end
        else
            frontCsv = fullfile(latestDir,'ga_front.csv');
            assert(exist(frontCsv,'file')==2, 'ga_front.csv bulunamadı.');
            Tbest = readtable(frontCsv);
            f1v = Tbest.f1; f2v = Tbest.f2;
            f1n = (f1v - min(f1v)) ./ max(eps, (max(f1v)-min(f1v)));
            f2n = (f2v - min(f2v)) ./ max(eps, (max(f2v)-min(f2v)));
            [~,row] = min(hypot(f1n, f2n));
        end
        gainsPF = struct('g_lo',Tbest.g_lo(row), 'g_mid',Tbest.g_mid(row), ...
                         'g_hi',Tbest.g_hi(row), 'PF_tau',Tbest.PF_tau(row), ...
                         'PF_gain',Tbest.PF_gain(row));
    end
end
if nargin < 4 || isempty(opts), opts = struct(); end

%% ------------------------------------------------------------------------
%% 2) VARSAYILANLARIN AYARLANMASI
getd = @(s,f,def) getfield_def(s,f,def); % kısa isim

dp_allow_frac = getd(opts,'dp_allow_frac',0.60);   % Δp_kv hedef / Δp_cap
alpha_lam     = getd(opts,'alpha_lam',0.15);       % laminer basınç payı
n_orf_set     = getd(opts,'n_orf_set',[5 6]);      % delik sayısı seçenekleri
dmm_step      = getd(opts,'dmm_step',0.05);        % d_o grid adımı [mm]
Cd_init       = getd(opts,'Cd_init', getd(getd(params,'orf',struct()),'CdInf',0.70));
mu_nom        = getd(opts,'mu_nom', getd(params,'mu_ref',0.9));

%% ------------------------------------------------------------------------
%% 3) KAZANÇLARIN VE PF’NİN SABİTLENMESİ
P = params;
nStories = size(P.M,1)-1;            % kat sayısı
% kazanc vektörü (alt 3 kat, üst 2 kat)
tg = ones(nStories,1) * gainsPF.g_mid;
loN = min(3,nStories); if loN>0, tg(1:loN) = gainsPF.g_lo; end
hiN = min(2,nStories); if hiN>0, tg(end-hiN+1:end) = gainsPF.g_hi; end
P.toggle_gain = tg;
% PF ayarları
if isfield(P,'cfg') && isfield(P.cfg,'PF')
    P.cfg.PF.tau  = gainsPF.PF_tau;
    P.cfg.PF.gain = gainsPF.PF_gain;
end

%% ------------------------------------------------------------------------
%% 4) SABİT KAZANÇLI SİMÜlASYON
O = struct('do_export',false,'quiet',true,'thermal_reset','each','order','natural', ...
           'use_orifice',true,'use_thermal',true, ...
           'mu_factors',[0.75 1.00 1.25], 'mu_weights',[0.2 0.6 0.2], 'thr', []);
S_worst = run_batch_windowed(scaled, P, O);

%% ------------------------------------------------------------------------
%% 5) Q95 VE BASINÇ LIMITLERİ
Q95_worst = NaN;
try
    vars = S_worst.table.Properties.VariableNames;
    if ismember('Q_q95_worst', vars)
        v = S_worst.table.Q_q95_worst; Q95_worst = max(v(:));
    elseif ismember('Q_q95_w', vars)
        v = S_worst.table.Q_q95_w;     Q95_worst = max(v(:));
    end
end

% model limiti veya thr üzerinden Δp_cap
dp_cap = getd(getd(O,'thr',struct()),'dP95_max',1.0e9);
dp_kv_target = dp_allow_frac * dp_cap;

%% ------------------------------------------------------------------------
%% 6) İLK ORİFİS BOYUTLANDIRMA (TÜRBÜLANS)
rho = P.rho;
Ao_req = abs(Q95_worst) / max(Cd_init,eps) * sqrt(rho / max(2*dp_kv_target,eps));
% d_o grid fonksiyonu (mm)
dgrid = @(dmm) (dmm_step * round(dmm./max(dmm_step,eps)));
best = struct('n_orf',NaN,'d_o',NaN,'A_o',NaN,'Cd',Cd_init,'Lori',NaN);
best_err = inf;
for n_orf = n_orf_set(:)'
    d_o  = sqrt(4*Ao_req/(pi*n_orf));       % m
    dmm  = dgrid(1000*d_o); d_o_q = dmm/1000; % kuantize m
    A_o  = n_orf * (pi*d_o_q^2/4);
    err  = abs(A_o - Ao_req)/max(A_o,eps);
    if err < best_err
        best_err = err; best.n_orf = n_orf; best.d_o = d_o_q; best.A_o = A_o;
    end
end

%% ------------------------------------------------------------------------
%% 7) LAMİNER KOL UZUNLUĞU
L_orif = alpha_lam * dp_kv_target * (pi * best.d_o^4) * best.n_orf / ...
         max(128*mu_nom*abs(Q95_worst),eps);
best.Lori = L_orif;

%% ------------------------------------------------------------------------
%% 8) Re VE Cd İTERASYONU
Cd = Cd_init;
for it = 1:2
    Q_hole = Q95_worst / max(best.n_orf,1);
    Re = 4*rho*abs(Q_hole) / max(pi*mu_nom*best.d_o,eps);
    Cd = P.orf.CdInf - (P.orf.CdInf - P.orf.Cd0) / ...
         (1 + (Re/max(P.orf.Rec,eps))^P.orf.p_exp);
    Ao_req = abs(Q95_worst) / max(Cd,eps) * sqrt(rho / max(2*dp_kv_target,eps));
    d_o = sqrt(4*Ao_req/(pi*best.n_orf));
    dmm = dgrid(1000*d_o); best.d_o = dmm/1000;
    best.A_o = best.n_orf*(pi*best.d_o^2/4);
end
best.Cd = Cd;

%% ------------------------------------------------------------------------
%% 9) GÜVENLİ DEBİ LİMİTİ
Qcap_big = 1.2 * abs(Q95_worst);

%% ------------------------------------------------------------------------
%% 10) PARAMETRE YAPISININ GÜNCELLENMESİ
P_sized = P;
P_sized.n_orf    = best.n_orf;
P_sized.orf.d_o  = best.d_o;
P_sized.A_o      = best.A_o;
P_sized.Qcap_big = Qcap_big;

%% ------------------------------------------------------------------------
%% 11) BOYUTLANDIRMA PAKETİ
sizing = struct();
sizing.Q95_worst    = Q95_worst;
sizing.dp_cap       = dp_cap;
sizing.dp_kv_target = dp_kv_target;
sizing.n_orf        = best.n_orf;
sizing.d_o          = best.d_o;
sizing.A_o          = best.A_o;
sizing.Cd           = best.Cd;
sizing.L_orif       = best.Lori;
sizing.Qcap_big     = Qcap_big;
sizing.alpha_lam    = alpha_lam;
sizing.notes        = 'Deterministik boyutlandırma; PF/kazanc sabit.';

%% ------------------------------------------------------------------------
%% 12) RAPOR VE CSV ÜRETİMİ
try
    [updates, T] = sizing_param_diff(params, P_sized, gainsPF, sizing, tg);
    fprintf('parametre.m için önerilen satırlar:\n');
    for i=1:numel(updates.lines)
        fprintf('  %s\n', updates.lines{i});
    end
    fprintf('  %% toggle_gain vektörü:\n  toggle_gain = %s;  %% (n-1)x1\n', mat2str(updates.toggle_gain_vec,3));
catch ME
    warning('sizing_param_diff başarısız: %s', ME.message);
end

end % make_sizing_case ana fonksiyon sonu

%% ========================================================================
%% ALT FONKSİYON: PARAMETRE FARKLARI
function [updates, T] = sizing_param_diff(P_old, P_sized, gainsPF, sizing, new_tg)
%SIZING_PARAM_DIFF Eski ve yeni parametreler arasındaki farkları hesaplar.
%   updates.lines parametre.m dosyasına yapıştırılacak satırları
%   içerir, T ise CSV tablo verisidir.

%% --- Eski değerleri çıkart ---
old_d_o   = getfield_def(P_old,'d_o',  getfield_def(getfield_def(P_old,'orf',struct()),'d_o',NaN));
old_Lori  = getfield_def(P_old,'Lori', getfield_def(getfield_def(P_old,'orf',struct()),'L_orif',NaN));
old_n_orf = getfield_def(P_old,'n_orf',NaN);
if ~isfinite(old_n_orf)
    old_Ao_tmp = getfield_def(P_old,'A_o',NaN);
    if isfinite(old_Ao_tmp) && isfinite(old_d_o)
        old_n_orf = old_Ao_tmp / (pi*(old_d_o^2)/4);
    end
end
if ~isfinite(old_Lori)
    c_lam0 = getfield_def(P_old,'c_lam0',NaN);
    mu_ref = getfield_def(P_old,'mu_ref',NaN);
    Ap     = getfield_def(P_old,'Ap',NaN);
    if all(isfinite([c_lam0 mu_ref Ap old_d_o]))
        old_Lori = c_lam0*(old_d_o^4)/(12*mu_ref*(Ap^2));
    end
end
if isfinite(old_n_orf) && isfinite(old_d_o)
    old_Ao = old_n_orf*pi*(old_d_o^2)/4;
else
    old_Ao = getfield_def(P_old,'A_o',NaN);
end
old_Qcap  = getfield_def(P_old,'Qcap_big',NaN);
old_tau   = getfield_def(getfield_def(getfield_def(P_old,'cfg',struct()),'PF',struct()),'tau',NaN);
old_gain  = getfield_def(getfield_def(getfield_def(P_old,'cfg',struct()),'PF',struct()),'gain',NaN);
old_tg    = getfield_def(P_old,'toggle_gain',[]);

%% --- Yeni değerleri çıkart ---
new_d_o   = getfield_def(getfield_def(P_sized,'orf',struct()),'d_o', getfield_def(P_sized,'d_o',NaN));
new_Lori  = getfield_def(getfield_def(P_sized,'orf',struct()),'L_orif', getfield_def(P_sized,'Lori', getfield_def(sizing,'L_orif',NaN)));
new_n_orf = getfield_def(P_sized,'n_orf',NaN);
if isfinite(new_n_orf) && isfinite(new_d_o)
    new_Ao = new_n_orf*pi*(new_d_o^2)/4;
else
    new_Ao = getfield_def(P_sized,'A_o',NaN);
end
new_Qcap  = getfield_def(P_sized,'Qcap_big',NaN);
new_tau   = getfield_def(gainsPF,'PF_tau',NaN);
new_gain  = getfield_def(gainsPF,'PF_gain',NaN);
new_tg    = new_tg(:);  % dışarıdan gelen yeni kazanc vektörü

%% --- toggle_gain eski değer ---
try
    nStories = size(P_old.M,1)-1;
catch
    nStories = numel(old_tg); if isempty(nStories) || nStories==0, nStories = numel(new_tg); end
end
old_tg = old_tg(:);
if isempty(old_tg)
    old_tg = ones(nStories,1)*getfield_def(gainsPF,'g_mid',1);
elseif isscalar(old_tg)
    old_tg = repmat(old_tg,nStories,1);
end

%% --- skaler eşitlik denetleyicisi ---
function tf = eqnum(a,b)
    if ~(isnumeric(a) && isnumeric(b)) || isempty(a) || isempty(b)
        tf = false; return;
    end
    a1 = a(1); b1 = b(1);
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
 'cfg.PF.gain (-)',   old_gain,  new_gain,  'cfg.PF.gain = %.6g;'
};

Name = {}; Old = {}; New = {}; Template = {};
for i=1:size(rows,1)
    if changed(rows{i,2}, rows{i,3})
        Name{end+1,1} = rows{i,1};
        if contains(rows{i,1},'n_orf')
            Old{end+1,1} = sprintf('%d', round(rows{i,2}));
            New{end+1,1} = sprintf('%d', round(rows{i,3}));
        else
            Old{end+1,1} = sprintf('%.6g', rows{i,2});
            New{end+1,1} = sprintf('%.6g', rows{i,3});
        end
        Template{end+1,1} = rows{i,4};
    end
end
T = table(Name, Old, New, Template);

%% --- toggle_gain farkı ---
tg_changed = numel(old_tg)~=numel(new_tg) || any(abs(old_tg(:)-new_tg(:))>1e-6);
if tg_changed
    T = [T; table({"toggle_gain ((n-1)x1)"}, {mat2str(old_tg,3)}, {mat2str(new_tg,3)}, {''}, 'VariableNames', T.Properties.VariableNames)];
    fprintf('\n[update] toggle_gain:\n  eski: %s\n  yeni: %s\n', mat2str(old_tg,3), mat2str(new_tg,3));
end

%% --- konsol çıktısı ---
fprintf('\n=== PARAMETRE GÜNCELLEME ÖNERİLERİ (sadece değişenler) ===\n');
for i=1:height(T)
    nm = string(T.Name{i});
    if T.Template{i}~=""
        if contains(nm,'n_orf')
            fprintf(' - %-18s : %s  -->  %s   (parametre.m: %s)\n', nm, T.Old{i}, T.New{i}, sprintf(T.Template{i}, round(str2double(T.New{i}))));
        else
            fprintf(' - %-18s : %s  -->  %s   (parametre.m: %s)\n', nm, T.Old{i}, T.New{i}, sprintf(T.Template{i}, str2double(T.New{i})));
        end
    else
        fprintf(' - %-18s : %s  -->  %s   (vektör)\n', nm, T.Old{i}, T.New{i});
    end
end
fprintf('===========================================================\n\n');

%% --- CSV Çıktısı ---
try
    outdir = fullfile('out'); if ~exist(outdir,'dir'), mkdir(outdir); end
    writetable(T, fullfile(outdir,'sizing_updates.csv'));
catch
end

%% --- parametre.m için satırlar ---
updates = struct(); updates.lines = {};
for i=1:height(T)
    if T.Template{i}~=""
        if contains(T.Name{i},'n_orf')
            updates.lines{end+1} = sprintf(T.Template{i}, round(str2double(T.New{i})));
        else
            updates.lines{end+1} = sprintf(T.Template{i}, str2double(T.New{i}));
        end
    end
end
updates.toggle_gain_vec = new_tg;
end

%% ========================================================================
%% YARDIMCI FONKSİYON: GÜVENLİ GETFIELD
function v = getfield_def(s, f, def)
%GETFIELD_DEF Yapı içinden alanı güvenli şekilde çeker, yoksa varsayılanı döner.
if nargin < 3, def = []; end
if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = def;
end
end

