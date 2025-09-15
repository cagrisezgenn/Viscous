function [X,F,gaout] = run_ga_driver(scaled, params, optsEval, optsGA)
%RUN_GA_DRIVER Hibrit GA sürücüsü
%   [X,F,GAOUT] = RUN_GA_DRIVER(SCALED, PARAMS, OPTSEVAL, OPTSGA)
% `scaled` ve `params` doğrudan argüman olarak verilebilir. Eğer boş
% bırakılırsa, gerekli veri seti ve parametreler `prepare_inputs` yardımcı
% fonksiyonu ile hazırlanır.

narginchk(0,4);

if nargin < 1 || isempty(scaled),  scaled  = []; end
if nargin < 2 || isempty(params),  params  = []; end
if nargin < 3 || isempty(optsEval), optsEval = struct; end
if nargin < 4 || isempty(optsGA),   optsGA   = struct; end

% Gerekli girdiler sağlanmadıysa otomatik hazırla
if isempty(scaled) || isempty(params)
    [scaled, params] = prepare_inputs(optsGA);
end
% === Parpool Açılışı (temizlik + iş parçacığı sınırı) ===
usePool = true;
try
    usePool = parpool_hard_reset(16);
catch ME
    warning('run_ga_driver:parpool', 'Parallel pool unavailable: %s', ME.message);
    usePool = false;
end
%RUN_GA_DRIVER Hibrit GA sürücüsü: önceden hazırlanmış `scaled` veri kümesi
% ve `params` yapısını kabul eder.
%   [X,F,GAOUT] = RUN_GA_DRIVER(SCALED, PARAMS, OPTSEVAL, OPTSGA)
% `scaled` veya `params` boş bırakılırsa `prepare_inputs` devreye girer ve
% gerekli verileri oluşturur. Uygunluk hesapları sırasında IO yapılmaz;
% yeniden ölçekleme yoktur.

meta = struct('thr', Utils.default_qc_thresholds(struct()));

assert(~isempty(scaled), 'run_ga_driver: scaled dataset is empty.');
assert(~isempty(params), 'run_ga_driver: params is empty.');

% ---------- Varsayılan değerlendirme ayarları (IO yok) ----------
optsEval.thr = Utils.default_qc_thresholds(Utils.getfield_default(optsEval,'thr', meta.thr));
%% GA Kurulumu
% GA amaç fonksiyonu ve optimizasyon seçeneklerini hazırla.
    rng(42);

    % Karar vektörü: [d_o_mm, n_orf, PF_tau, PF_gain, Cd0, CdInf, p_exp, Lori_mm, hA_W_perK, Dp_mm, d_w_mm, D_m_mm, n_turn, mu_ref]
    lb = [1.0, 8, 0.4, 2.5, 0.60, 0.75, 0.90, 120, 600, 120, 10,  90,  6, 0.80, 0.0];
ub = [3.0,8, 0.90, 5, 0.90, 1.00, 1.50, 200, 600, 240, 16, 160, 18, 2.00, 3];
    IntCon = [2 13];  % n_orf ve n_turn tam sayı

    % Veri seti imzası üret (önbellek anahtarı)
        dsig = sum([scaled.IM]) + sum([scaled.PGA]);
    optsEval.dsig = dsig;
    obj = @(x) eval_design_fast(x, scaled, params, optsEval); % içerde kuantize/clamplar

    options = optimoptions('gamultiobj', ...
       'PopulationSize',    Utils.getfield_default(optsGA,'PopulationSize',40), ...
       'MaxGenerations',    Utils.getfield_default(optsGA,'MaxGenerations',4), ...
       'CrossoverFraction', Utils.getfield_default(optsGA,'CrossoverFraction',0.9), ...
       'MutationFcn',       Utils.getfield_default(optsGA,'MutationFcn',{@mutationgaussian,0.2,0.4}), ...
       'ParetoFraction',    Utils.getfield_default(optsGA,'ParetoFraction',0.7), ...
       'StallGenLimit',     Utils.getfield_default(optsGA,'StallGenLimit',150), ...
       'DistanceMeasureFcn','distancecrowding', ...
       'OutputFcn',         @(options,state,flag) ga_out_best_pen(options,state,flag, scaled, params, optsEval), ...
       'UseParallel',       usePool && Utils.getfield_default(optsGA,'UseParallel',true), ...
       'Display','iter','PlotFcn',[], 'FunctionTolerance',1e-5);

    %% Başlangıç Popülasyonu
    % Izgaraya hizalı ilk popülasyonu oluştur (tohumlarla birlikte).
        if isempty(Utils.getfield_default(optsGA,'InitialPopulationMatrix',[]))
            step_vec = [0.1 NaN 0.01 0.02 0.01 0.01 0.05 1 25 1 0.5 5 NaN 0.05 0.1];
            P0 = Utils.initial_pop_grid(lb, ub, options.PopulationSize, step_vec);
            % Tohumlar (yeni sınırlar içinde uygulanabilir)
           seed = [
 2.50  4  0.03 0.60 0.50 0.75 0.80  60 200  95  7  55  3 0.40;   % alt sınırların köşesi
 3.90  8  0.12 1.90 0.75 1.00 1.50 160 900 195 18 125 12 1.50;   % üst sınırların köşesi
 3.20  6  0.08 1.25 0.60 0.85 1.15 100 550 145 12  90  8 1.00;   % ortalama değerler
 2.80  4  0.04 0.90 0.55 0.80 1.00  80 300 120  9  70  6 0.70;   % düşük n_orf / n_turn
 3.40  7  0.09 1.60 0.65 0.90 1.30 140 700 170 16 110 10 1.30;   % yüksek PF_gain ve n_turn
 2.60  5  0.06 1.10 0.60 0.82 1.20  50 250 100 10  60  5 0.90;   % düşük Lori_mm, orta diğerleri
];
            % Phase 7: if decision dim expanded (PF_t_on), pad seeds
            if size(seed,2) < size(P0,2)
                seed = [seed, 0.75*ones(size(seed,1), size(P0,2)-size(seed,2))];
            end
            % Tohumları güvenli aralığa sıkıştır (lb/ub)
            seed = max(min(seed, ub), lb);
            ns = min(size(seed,1), size(P0,1));
            P0(1:ns,:) = seed(1:ns,:);
            options = optimoptions(options,'InitialPopulationMatrix', P0);
        end

    [X,F,exitflag,output,population,scores] = gamultiobj(obj, numel(lb), [],[],[],[], lb, ub, [], IntCon, options);
    gaout = struct('exitflag',exitflag,'output',output);

    %% Sonuçların Paketlenmesi
    % GA tamamlandıktan sonra sonuçları dosyalara kaydet.
    tstamp = datestr(now,'yyyymmdd_HHMMSS_FFF');
    outdir = fullfile('out', ['ga_' tstamp]);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    date_str = tstamp;
    front = struct('X',X,'F',F,'options',options,'meta',meta,'date_str',date_str);
    save(fullfile(outdir,'ga_front.mat'), '-struct', 'front', '-v7.3');

    % === Re-evaluate Pareto designs to collect metrics per-row ===
    nF = size(X,1);
    x10_max_damperli  = zeros(nF,1);  a10abs_max_damperli = zeros(nF,1);
    dP95   = zeros(nF,1);  Qcap95 = zeros(nF,1); cav_pct = zeros(nF,1);
    T_end  = zeros(nF,1);  mu_end  = zeros(nF,1);
    PF_p95 = zeros(nF,1);  Q_q50 = zeros(nF,1);  Q_q95 = zeros(nF,1);
    dP50   = zeros(nF,1);  Toil = nan(nF,1); Tsteel = nan(nF,1);
    energy_tot_sum   = zeros(nF,1);  E_orifice_sum = zeros(nF,1);   E_struct_sum  = zeros(nF,1); E_ratio = zeros(nF,1); P_mech_sum = zeros(nF,1);
    PFA_mean   = zeros(nF,1);  IDR_mean  = zeros(nF,1);

    % ceza bileşenleri (eval ile aynı)
    pen     = zeros(nF,1);
    pen_dP  = zeros(nF,1); pen_Qcap = zeros(nF,1); pen_cav = zeros(nF,1); pen_T = zeros(nF,1); pen_mu = zeros(nF,1);
    lambda  = Utils.getfield_default(optsEval,'penalty_scale',5);   % daha yumuşak ceza
    pwr     = Utils.getfield_default(optsEval,'penalty_power',1.0);
    W       = Utils.getfield_default(optsEval,'penalty_weights', struct('dP',1,'Qcap',1,'cav',1.5,'T',0.5,'mu',0.5));
    rel = @(v,lim) max(0,(v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0,(lim - v)./max(lim,eps)).^pwr;

    Opost = struct('thr', meta.thr);

    parfor i = 1:nF
        tg = @(tbl,var,def) Utils.table_get(tbl,var,def);
        Xi = quant_clamp_x(X(i,:));
        Pi = decode_params_from_x(params, Xi);
        Si = run_batch_windowed(scaled, Pi, Opost);

        % Tablo sütunlarını güvenle al
            v_PFA = Si.table.PFA;
            v_IDR = Si.table.IDR;
            v_x10 = Si.table.x10_max_D;
            v_a10 = Si.table.a10abs_max_D;

            v_dP95 = Si.table.dP95;
            v_Qcap = Si.table.Qcap95;
            v_cav  = Si.table.cav_pct;
            v_T_end = Si.table.T_end;
            v_mu   = Si.table.mu_end;

            v_PF_p95 = tg(Si.table,'PF_p95',0);
            v_Q_q50  = tg(Si.table,'Q_q50',0);
            v_Q_q95  = tg(Si.table,'Q_q95',0);
            v_dP50 = tg(Si.table,'dP50',0);

        v_E_orifice_sum = tg(Si.table,'E_orifice_sum',0);
        v_E_struct_sum = tg(Si.table,'E_struct_sum',0);
        v_P_mech_sum  = tg(Si.table,'P_mech_sum',0);

        PFA_mean(i) = mean(v_PFA(:));
        IDR_mean(i) = mean(v_IDR(:));

        % Aggregate across records (dataset) to scalars per design
        x10_max_damperli(i)  = max(v_x10(:));
        a10abs_max_damperli(i)  = max(v_a10(:));

        dP95(i)   = max(v_dP95(:));
        Qcap95(i) = max(v_Qcap(:));
        cav_pct(i)   = max(v_cav(:));
        T_end(i)   = max(v_T_end(:));
        if isempty(v_mu), v_mu = 1; end
        mu_end(i)  = min(v_mu(:));

        PF_p95(i)  = max(v_PF_p95(:));
        Q_q50(i)   = max(v_Q_q50(:));
        Q_q95(i)   = max(v_Q_q95(:));
        dP50(i)  = max(v_dP50(:));
        Toil(i)   = max(tg(Si.table,'T_oil_end',NaN));
        Tsteel(i) = max(tg(Si.table,'T_steel_end',NaN));

        E_orifice_sum(i)    = sum(v_E_orifice_sum(:));
        E_struct_sum(i)   = sum(v_E_struct_sum(:));
        energy_tot_sum(i)   = E_orifice_sum(i) + E_struct_sum(i);
        E_ratio(i) = (E_struct_sum(i)>0) * (E_orifice_sum(i)/max(E_struct_sum(i),eps));
        P_mech_sum(i)  = sum(v_P_mech_sum(:));

        % Penalty parts (mirror eval_design_fast): mean over records
        vv_dP   = v_dP95(:); if isempty(vv_dP), vv_dP = 0; end
        vv_Qcap = v_Qcap(:); if isempty(vv_Qcap), vv_Qcap = 0; end
        vv_cav  = v_cav(:);  if isempty(vv_cav), vv_cav = 0; end
        vv_T_end = v_T_end(:); if isempty(vv_T_end), vv_T_end = 0; end
        vv_mu   = v_mu(:);   if isempty(vv_mu), vv_mu = 1; end

        pen_dP(i)   = mean(rel(vv_dP,   Opost.thr.dP95_max));
        pen_Qcap(i) = mean(rel(vv_Qcap, Opost.thr.Qcap95_max));
        if Opost.thr.cav_pct_max<=0
            pen_cav(i) = mean(max(0,vv_cav).^pwr);
        else
            pen_cav(i) = mean(rel(vv_cav, Opost.thr.cav_pct_max));
        end
        pen_T(i)    = mean(rel(vv_T_end, Opost.thr.T_end_max));
        pen_mu(i)   = mean(rev(vv_mu,   Opost.thr.mu_end_min));
        pen(i)      = lambda*(W.dP*pen_dP(i)+W.Qcap*pen_Qcap(i)+W.cav*pen_cav(i)+W.T*pen_T(i)+W.mu*pen_mu(i));
    end

    % Satır başına dizilerden T tablosunu oluştur
    data = [X F PFA_mean IDR_mean pen pen_dP pen_Qcap pen_cav pen_T pen_mu ...
            x10_max_damperli a10abs_max_damperli dP95 Qcap95 cav_pct T_end mu_end PF_p95 ...
            Q_q50 Q_q95 dP50 Toil Tsteel energy_tot_sum E_orifice_sum E_struct_sum E_ratio P_mech_sum];
    T = array2table(data, 'VariableNames', ...
       {'d_o_mm','n_orf','PF_tau','PF_gain','Cd0','CdInf','p_exp','Lori_mm','hA_W_perK','Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref','PF_t_on', ...
        'f1','f2','PFA_mean','IDR_mean','pen','pen_dP','pen_Qcap','pen_cav','pen_T','pen_mu', ...
        'x10_max_damperli','a10abs_max_damperli','dP95','Qcap95','cav_pct','T_end','mu_end','PF_p95', ...
        'Q_q50','Q_q95','dP50','T_oil_end','T_steel_end','energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio','P_mech_sum'});

    % === BASELINE (pre-GA) ROW: params başlangıcıyla tek koşu, ilk satır ===
    T = prepend_baseline_row(T, params, scaled, Opost, lambda, pwr, W);

    write_pareto_results(T, outdir);

    % En iyi K tasarımın parametre listesi
    K = min(10, size(X,1));
    if K > 0
        idx = unique(round(linspace(1,size(X,1),K)));
        params_list = cell(numel(idx),1);
        for i = 1:numel(idx)
            params_list{i} = decode_params_from_x(params, X(idx(i),:)); %#ok<AGROW>
        end
        tmp_struct = struct('params_list',{params_list});
        save(fullfile(outdir,'ga_front.mat'), '-struct', 'tmp_struct', '-append');
    end

end

function [f, meta] = eval_design_fast(x, scaled, params_base, optsEval)
    % ızgaralara oturt
    x = x(:)';
    % kuantize et
    x(1) = Utils.quantize_step(x(1),0.05);  % d_o_mm
    x(3) = Utils.quantize_step(x(3),0.01);  % PF_tau
    x(4) = Utils.quantize_step(x(4),0.02);  % PF_gain
    x(5) = Utils.quantize_step(x(5),0.01);  % Cd0
    x(6) = Utils.quantize_step(x(6),0.01);  % CdInf
    x(7) = Utils.quantize_step(x(7),0.05);  % p_exp
    if numel(x) >= 8,  x(8)  = Utils.quantize_step(x(8),1);   end % Lori_mm
    if numel(x) >= 9,  x(9)  = Utils.quantize_step(x(9),25);  end % hA_W_perK
    if numel(x) >=10,  x(10) = Utils.quantize_step(x(10),1);  end % Dp_mm
    if numel(x) >=11,  x(11) = Utils.quantize_step(x(11),0.5);end % d_w_mm
    if numel(x) >=12,  x(12) = Utils.quantize_step(x(12),5);  end % D_m_mm
    if numel(x) >=13,  x(13) = round(x(13));                end % n_turn
    if numel(x) >=14,  x(14) = Utils.quantize_step(x(14),0.05); end
    if numel(x) >=15,  x(15) = Utils.quantize_step(x(15),0.25); end % mu_ref
    x(2)  = round(max(x(2),1));        % n_orf tam sayı, >=1
    if numel(x) >=13, x(13) = round(max(x(13),1)); end

    

    persistent memo;
    if isempty(memo), memo = containers.Map(); end
    % Bellek anahtarlarının farklı veri setleri arasında çakışmaması için tuz ekle
    dsig = 0;
        dsig = sum([scaled.IM]) + sum([scaled.PGA]);
    key = jsonencode([x, dsig]);
    if isKey(memo, key)
        meta = memo(key); f = meta.f; return;
    end

    P = decode_params_from_x(params_base, x);

    O = struct();
    if nargin >= 4 && ~isempty(optsEval), O = optsEval; end
    % policy/order vars referenced by run_batch_windowed default to each/natural

    % Güvenli değerlendirme (GA sırasında IO yok)
        S = run_batch_windowed(scaled, P, O);

    % --- HARD FILTER (erken eleme) ---
    dP95v = S.table.dP95;
    qcapv = S.table.Qcap95;
    cavv  = S.table.cav_pct;
    if any(dP95v > 1e9) || any(qcapv > 0.90) || any(cavv > 0.01)
        f = [1e6, 1e6];
        meta = struct('x',x,'f',f,'hard_kill',true);
        % --- penalties: ensure numeric zeros for hard-kill
        meta.pen      = 0;
        meta.pen_dP   = 0;
        meta.pen_Qcap = 0;
        meta.pen_cav  = 0;
        meta.pen_T    = 0;
        meta.pen_mu   = 0;
        memo_store('set', jsonencode([x, dsig]), meta);
        return;
    end

    f1 = mean(S.table.PFA);
    f2 = mean(S.table.IDR);
    thr = O.thr;
    dP95v   = S.table.dP95;
    qcapv   = S.table.Qcap95;
    cavv    = S.table.cav_pct;
    T_endv   = S.table.T_end;
    mu_endv  = S.table.mu_end;
    % --- PENALTY (soft, multiplicative) ---
    lambda  = Utils.getfield_default(optsEval,'penalty_scale',10);
    pwr     = Utils.getfield_default(optsEval,'penalty_power',1.0);
    rel = @(v,lim) max(0, (v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0, (lim - v)./max(lim,eps)).^pwr;
    W = Utils.getfield_default(optsEval,'penalty_weights', struct('dP',1,'Qcap',1,'cav',2,'T',0.5,'mu',0.5));
    pen_dP   = isfield(thr,'dP95_max')   * mean(rel(dP95v,  thr.dP95_max));
    pen_Qcap = isfield(thr,'Qcap95_max') * mean(rel(qcapv,  thr.Qcap95_max));
    if isfield(thr,'cav_pct_max')
        if thr.cav_pct_max<=0, pen_cav = mean(max(0,cavv).^pwr); else, pen_cav = mean(rel(cavv, thr.cav_pct_max)); end
    else
        pen_cav = 0;
    end
    pen_T    = isfield(thr,'T_end_max')  * mean(rel(T_endv,  thr.T_end_max));
    pen_mu   = isfield(thr,'mu_end_min') * mean(rev(mu_endv, thr.mu_end_min));
    pen = W.dP*pen_dP + W.Qcap*pen_Qcap + W.cav*pen_cav + W.T*pen_T + W.mu*pen_mu;
    % multiplicative penalty keeps scale of objectives
    f = [f1, f2] .* (1 + lambda*pen);
    % Build meta and enforce numeric penalties
    meta = struct('x',x,'f',f,'PFA_mean',f1,'IDR_mean',f2);
    Penalty   = lambda*pen;
    pen_parts = struct('dP',pen_dP,'Qcap',pen_Qcap,'cav',pen_cav,'T',pen_T,'mu',pen_mu);
    % --- penalties: force numeric (no NaNs)
    if ~exist('Penalty','var') || ~isfinite(Penalty), Penalty = 0; end
    if ~exist('pen_parts','var') || ~isstruct(pen_parts), pen_parts = struct(); end
    pf = {'dP','Qcap','cav','T','mu'};
    for ii=1:numel(pf)
        fn = pf{ii};
        if ~isfield(pen_parts,fn) || ~isfinite(pen_parts.(fn)), pen_parts.(fn) = 0; end
    end
    meta.pen      = Penalty;
    meta.pen_dP   = pen_parts.dP;
    meta.pen_Qcap = pen_parts.Qcap;
    meta.pen_cav  = pen_parts.cav;
    meta.pen_T    = pen_parts.T;
    meta.pen_mu   = pen_parts.mu;
    meta.pen_parts = pen_parts;
    % --- append damperli peak metrics (zeros if unavailable)
    x10_max_damperli_local = 0; a10abs_max_damperli_local = 0;
            if isfield(S,'table') && istable(S.table)
                if ismember('x10_max_D', S.table.Properties.VariableNames)
                        x10_max_damperli_local = max(S.table.x10_max_D);
                end
                if ismember('a10abs_max_D', S.table.Properties.VariableNames)
                        a10abs_max_damperli_local = max(S.table.a10abs_max_D);
                end
            end
    meta.x10_max_damperli    = x10_max_damperli_local;
    meta.a10abs_max_damperli = a10abs_max_damperli_local;
    % === Penaltı sürücüleri ve diagnostikleri ekle (varsa) ===
        if isfield(S,'table') && istable(S.table)
            candCols = { ...
              'dP95','Qcap95','cav_pct','T_end','mu_end', ...
              'PF_p95', ...
              'Q_q50','Q_q95','dP50', ...
              'T_oil_end','T_steel_end', ...
              'energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio','P_mech_sum', ...
              'x10_max_D','a10abs_max_D','P_mech','E_orifice','E_struct','Re_max' ...
            };
            % Gerekirse alan adları için takma adlar sağla
                if ismember('E_orifice_sum', S.table.Properties.VariableNames) && ...
                   ismember('E_struct_sum', S.table.Properties.VariableNames) && ...
                   ~ismember('energy_tot_sum', S.table.Properties.VariableNames)
                    S.table.energy_tot_sum = S.table.E_orifice_sum + S.table.E_struct_sum; %#ok<AGROW>
                end
            for jj = 1:numel(candCols)
                nm = candCols{jj};
                if ismember(nm, S.table.Properties.VariableNames)
                    val = S.table.(nm);
                        if ~all(isfinite(val(:)))
                            val(~isfinite(val)) = 0;
                        end
                    meta.(nm) = val;
                end
            end
        end

memo(key) = meta;
memo_store('set', key, meta);

    end

function T = prepend_baseline_row(T, params, scaled, Opost, lambda, pwr, W)
    rel = @(v,lim) max(0,(v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0,(lim - v)./max(lim,eps)).^pwr;
        % Build X0 from params (with sensible fallbacks)
        X0 = nan(1,13);
        % d_o_mm and n_orf
            if isfield(params,'orf') && isfield(params.orf,'d_o') && ~isempty(params.orf.d_o)
                X0(1) = 1e3 * params.orf.d_o;
            else
                X0(1) = 3.00;
            end
            if isfield(params,'n_orf') && ~isempty(params.n_orf)
                X0(2) = params.n_orf;
            else
                X0(2) = 5;
            end
        % PF parameters
            if isfield(params,'cfg') && isfield(params.cfg,'PF')
                X0(3) = Utils.getfield_default(params.cfg.PF,'tau',1.0);
                X0(4) = Utils.getfield_default(params.cfg.PF,'gain',0.85);
            else
                X0(3:4) = [1.0 0.85];
            end

        % Orifis katsayıları
            X0(5) = Utils.getfield_default(params.orf,'Cd0',0.61);
            X0(6) = Utils.getfield_default(params.orf,'CdInf',0.80);
            X0(7) = Utils.getfield_default(params.orf,'p_exp',1.10);

        % Lori, hA ve geometri
            X0(8) = Utils.getfield_default(params,'Lori',0.10)*1e3;
            X0(9) = Utils.getfield_default(params.thermal,'hA_W_perK',450);
            X0(10) = Utils.getfield_default(params,'Dp',0.125)*1e3;
            X0(11) = Utils.getfield_default(params,'d_w',0.012)*1e3;
            X0(12) = Utils.getfield_default(params,'D_m',0.080)*1e3;
            X0(13) = Utils.getfield_default(params,'n_turn',8);
            X0(14) = Utils.getfield_default(params,'mu_ref',0.9);

        % Simulate baseline with same post-eval options
        X0 = quant_clamp_x(X0);
        P0    = decode_params_from_x(params, X0);
        S0    = run_batch_windowed(scaled, P0, Opost);
        T0bl  = S0.table;
        % Objectives
        f0 = [NaN NaN];
        [f0, ~] = eval_design_fast(X0, scaled, params, Opost);
            PFA0 = mean(T0bl.PFA);
            IDR0 = mean(T0bl.IDR);
        % Penalty parts (same formula)
        dP95_0     = max(T0bl.dP95);
        Qcap95_0   = max(T0bl.Qcap95);
        cav_pct_0  = max(T0bl.cav_pct);
        T_end_0    = max(T0bl.T_end);
        mu_end_0   = min(T0bl.mu_end);
        pen_dP_0   = rel(dP95_0,  Opost.thr.dP95_max);
        pen_Qcap_0 = rel(Qcap95_0,Opost.thr.Qcap95_max);
        if Opost.thr.cav_pct_max<=0, pen_cav_0 = max(0,cav_pct_0).^pwr; else, pen_cav_0 = rel(cav_pct_0,Opost.thr.cav_pct_max); end
        pen_T_0    = rel(T_end_0,  Opost.thr.T_end_max);
        pen_mu_0   = rev(mu_end_0, Opost.thr.mu_end_min);
        pen_0      = lambda*(W.dP*pen_dP_0 + W.Qcap*pen_Qcap_0 + W.cav*pen_cav_0 + W.T*pen_T_0 + W.mu*pen_mu_0);

        % Baz satırı için tablonun ilk satırını kopyala
        vn = T.Properties.VariableNames;
        T0 = T(1,:);
        T0{1,:} = nan(1,width(T0));

        % T0'ya değer yazarken sütun tiplerini koru

        % karar değişkenleri
        assign('d_o_mm',  X0(1));
        assign('n_orf',   X0(2));
        assign('PF_tau',  X0(3));
        assign('PF_gain', X0(4));
        assign('Cd0',     X0(5));
        assign('CdInf',   X0(6));
        assign('p_exp',   X0(7));
        assign('Lori_mm', X0(8));
        assign('hA_W_perK', X0(9));
        assign('Dp_mm',   X0(10));
        assign('d_w_mm',  X0(11));
        assign('D_m_mm',  X0(12));
        assign('n_turn',  X0(13));
        assign('mu_ref',  X0(14));

        % amaç fonksiyonları
        assign('f1', f0(1));
        assign('f2', f0(2));
        assign('PFA_mean', PFA0);
        assign('IDR_mean', IDR0);

        % cezalar
        assign('pen',     pen_0);
        assign('pen_dP',  pen_dP_0);
        assign('pen_Qcap',pen_Qcap_0);
        assign('pen_cav', pen_cav_0);
        assign('pen_T',   pen_T_0);
        assign('pen_mu',  pen_mu_0);

        % damper tepe değerleri
            if ismember('x10_max_damperli', vn) && ismember('x10_max_D', T0bl.Properties.VariableNames)
                assign('x10_max_damperli', max(T0bl.x10_max_D));
            end
            if ismember('a10abs_max_damperli', vn) && ismember('a10abs_max_D', T0bl.Properties.VariableNames)
                assign('a10abs_max_damperli', max(T0bl.a10abs_max_D));
            end

        % diğer diagnostikler
            if ismember('dP95', vn) && ismember('dP95', T0bl.Properties.VariableNames)
                assign('dP95', max(T0bl.dP95));
            end
            if ismember('Qcap95', vn) && ismember('Qcap95', T0bl.Properties.VariableNames)
                assign('Qcap95', max(T0bl.Qcap95));
            end
            if ismember('cav_pct', vn) && ismember('cav_pct', T0bl.Properties.VariableNames)
                assign('cav_pct', max(T0bl.cav_pct));
            end
            if ismember('T_end', vn) && ismember('T_end', T0bl.Properties.VariableNames)
                assign('T_end', max(T0bl.T_end));
            end
            if ismember('mu_end', vn) && ismember('mu_end', T0bl.Properties.VariableNames)
                assign('mu_end', min(T0bl.mu_end));
            end
            if ismember('PF_p95', vn) && ismember('PF_p95', T0bl.Properties.VariableNames)
                assign('PF_p95', max(T0bl.PF_p95));
            end
            if ismember('Q_q50', vn) && ismember('Q_q50', T0bl.Properties.VariableNames)
                assign('Q_q50', max(T0bl.Q_q50));
            end
            if ismember('Q_q95', vn) && ismember('Q_q95', T0bl.Properties.VariableNames)
                assign('Q_q95', max(T0bl.Q_q95));
            end
            if ismember('dP50', vn) && ismember('dP50', T0bl.Properties.VariableNames)
                assign('dP50', max(T0bl.dP50));
            end
            if ismember('T_oil_end', vn) && ismember('T_oil_end', T0bl.Properties.VariableNames)
                assign('T_oil_end', max(T0bl.T_oil_end));
            end
            if ismember('T_steel_end', vn) && ismember('T_steel_end', T0bl.Properties.VariableNames)
                if iscell(T.T_steel_end)
                    T0.T_steel_end = {max(T0bl.T_steel_end)};
                else
                    T0.T_steel_end = max(T0bl.T_steel_end);
                end
            end
            if ismember('E_orifice_sum', vn) && ismember('E_orifice_sum', T0bl.Properties.VariableNames)
                assign('E_orifice_sum', sum(T0bl.E_orifice_sum));
            end
            if ismember('E_struct_sum', vn) && ismember('E_struct_sum', T0bl.Properties.VariableNames)
                assign('E_struct_sum', sum(T0bl.E_struct_sum));
            end
            if ismember('energy_tot_sum', vn)
                if ismember('energy_tot_sum', T0bl.Properties.VariableNames)
                    assign('energy_tot_sum', sum(T0bl.energy_tot_sum));
                else
                        assign('energy_tot_sum', T0.E_orifice_sum + T0.E_struct_sum);
                end
            end
            if ismember('E_ratio', vn)
                    assign('E_ratio', (T0.E_struct_sum>0) * (T0.E_orifice_sum / max(T0.E_struct_sum, eps)));
            end
            if ismember('P_mech_sum', vn) && ismember('P_mech_sum', T0bl.Properties.VariableNames)
                assign('P_mech_sum', sum(T0bl.P_mech_sum));
            end

        % Bas satırı en üste ekle
        T = [T0; T];
    function assign(name, val)
        if ismember(name, vn)
            % Coerce val to numeric when appropriate
            if iscell(val)
                if numel(val)==1 && isnumeric(val{1})
                    val = val{1};
                elseif numel(val)==1 && ischar(val{1})
                    vtmp = str2double(val{1});
                    if isfinite(vtmp), val = vtmp; else, val = NaN; end
                else
                    % multi-cell: if destination is cell, pass through; else NaN
                    if ~iscell(T.(name))
                        val = NaN;
                    end
                end
            elseif ischar(val)
                vtmp = str2double(val);
                if isfinite(vtmp), val = vtmp; else, val = NaN; end
            end
            if iscell(T.(name))
                T0.(name) = {val};
            else
                if isnumeric(val)
                    T0.(name) = val;
                else
                    T0.(name) = NaN;
                end
            end
        end
    end
end

function write_pareto_results(T, outdir)
    writetable(T, fullfile(outdir,'ga_front.csv'));
end


function xq = quant_clamp_x(x)
    % Apply the same quantization/clamps as in eval
    xq = x;
    if isvector(x), xq = x(:)'; end
    xq(:,1) = Utils.quantize_step(xq(:,1),0.05);
    xq(:,3) = Utils.quantize_step(xq(:,3),0.01);
    xq(:,4) = Utils.quantize_step(xq(:,4),0.02);
    xq(:,5) = Utils.quantize_step(xq(:,5),0.01);
    xq(:,6) = Utils.quantize_step(xq(:,6),0.01);
    xq(:,7) = Utils.quantize_step(xq(:,7),0.05);
    if size(xq,2) >= 8,  xq(:,8)  = Utils.quantize_step(xq(:,8),1);   end
    if size(xq,2) >= 9,  xq(:,9)  = Utils.quantize_step(xq(:,9),25);  end
    if size(xq,2) >=10,  xq(:,10) = Utils.quantize_step(xq(:,10),1);  end
    if size(xq,2) >=11,  xq(:,11) = Utils.quantize_step(xq(:,11),0.5);end
    if size(xq,2) >=12,  xq(:,12) = Utils.quantize_step(xq(:,12),5);  end
    if size(xq,2) >=13,  xq(:,13) = round(xq(:,13));                end
    if size(xq,2) >=14,  xq(:,14) = Utils.quantize_step(xq(:,14),0.05); end
    if size(xq,2) >=15,  xq(:,15) = Utils.quantize_step(xq(:,15),0.25); end
    xq(:,2) = round(max(xq(:,2),1));
    if size(xq,2) >=13, xq(:,13) = round(max(xq(:,13),1)); end
    if isvector(x), xq = xq(:)'; end
end

function out = memo_store(cmd, key, val)
    persistent MEM;
    if isempty(MEM), MEM = containers.Map(); end
    switch lower(cmd)
        case 'get'
            if isKey(MEM,key), out = MEM(key); else, out = []; end
        case 'set'
            MEM(key) = val; out = true;
        case 'clear'
            MEM = containers.Map(); out = true;
        otherwise
            out = [];
    end
end

function P = decode_params_from_x(params_base_, x_)
    d_o_mm = x_(1); n_orf = round(x_(2));
    PF_tau = x_(3); PF_gain = x_(4);
    Cd0 = x_(5); CdInf = x_(6); p_exp = x_(7);
    Lori_mm = x_(8); hA_W_perK = x_(9);
    Dp_mm = x_(10); d_w_mm = x_(11); D_m_mm = x_(12);
    n_turn = round(x_(13)); mu_ref = x_(14);
    if numel(x_) >= 15, PF_t_on = x_(15); else, PF_t_on = NaN; end
    P = params_base_;
    P.orf.d_o = d_o_mm * 1e-3;         % mm'den m'ye
    P.n_orf   = n_orf;
    P.orf.Cd0   = Cd0;
    P.orf.CdInf = CdInf;
    P.orf.p_exp = p_exp;
    P.Lori   = Lori_mm * 1e-3;
    P.mu_ref = mu_ref;
    P.Dp     = Dp_mm * 1e-3;
    P.d_w    = d_w_mm * 1e-3;
    P.D_m    = D_m_mm * 1e-3;
    P.n_turn = n_turn;
    if isfield(P,'thermal')
        P.thermal.hA_W_perK = hA_W_perK;
    end
    P = build_params(P);
    if isfield(P,'cfg') && isfield(P.cfg,'PF')
        P.cfg.PF.tau  = PF_tau;
        P.cfg.PF.gain = PF_gain;
        if isfinite(PF_t_on)
            P.cfg.PF.t_on = PF_t_on;
            if isfield(P.cfg.PF,'auto_t_on') && P.cfg.PF.auto_t_on
                % GA vektörü t_on sağlıyorsa otomatik t_on'u devre dışı bırak
                P.cfg.PF.auto_t_on = false;
            end
        end
    end
end
function [state, options, optchanged] = ga_out_best_pen(options, state, flag, scaled, params, optsEval)
    optchanged = false;
    if strcmp(flag,'iter')
        [~, idx] = min(sum(state.Score,2));
        bestx = state.Population(idx,:);
        [f_curr, meta_curr] = eval_design_fast(bestx, scaled, params, optsEval);
        fprintf('Gen %d: f1=%g f2=%g pen=%g\n', state.Generation, f_curr(1), f_curr(2), meta_curr.pen);
    end
end
