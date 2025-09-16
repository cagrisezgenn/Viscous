function [X,F,gaout] = run_ga_driver(scaled, params, optsEval, optsGA)
%RUN_GA_DRIVER Hibrit GA sürücüsü
%   [X,F,GAOUT] = RUN_GA_DRIVER(SCALED, PARAMS, OPTSEVAL, OPTSGA)
% `scaled` ve `params` doğrudan argüman olarak verilebilir. Eğer boş
% bırakılırsa, gerekli veri seti ve parametreler `prepare_inputs`
% fonksiyonu ile hazırlanır. Uygunluk hesapları sırasında IO yapılmaz;
% yeniden ölçekleme beklenmez.
%   Bu sürüm, daha önce ayrı dosyalarda yer alan tüm yardımcı
%   fonksiyonları (hazırlık, kayıt döngüsü, dinamik çözüm vb.) aynı dosya
%   içinde toplayarak tek başına çalışabilir bir GA iş akışı sunar.

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

assert(~isempty(scaled), 'run_ga_driver: scaled dataset is empty.');
assert(~isempty(params), 'run_ga_driver: params is empty.');

% ---------- Varsayılan değerlendirme ayarları (IO yok) ----------
optsEval.thr = util_default_qc_thresholds(util_getfield_default(optsEval,'thr', struct()));
%% GA Kurulumu
% GA amaç fonksiyonu ve optimizasyon seçeneklerini hazırla.
    rng(42);

    % Karar vektörü: [d_o_mm, n_orf, PF_tau, PF_gain, Cd0, CdInf, p_exp, Lori_mm, hA_W_perK, Dp_mm, d_w_mm, D_m_mm, n_turn, mu_ref, PF_t_on]
    lb = [1.0, 4, 0.40, 2.5, 0.60, 0.75, 0.90, 120, 600, 120, 10,  90,  6, 0.80, 0.0];
    ub = [3.0, 8, 0.90, 5.0, 0.90, 1.00, 1.50, 200, 600, 240, 16, 160, 18, 2.00, 3.0];
    IntCon = [2 13];  % n_orf ve n_turn tam sayı

    obj = @(x) eval_design_fast(x, scaled, params, optsEval); % içerde kuantize/clamplar

    options = optimoptions('gamultiobj', ...
       'PopulationSize',    util_getfield_default(optsGA,'PopulationSize',40), ...
       'MaxGenerations',    util_getfield_default(optsGA,'MaxGenerations',4), ...
       'CrossoverFraction', util_getfield_default(optsGA,'CrossoverFraction',0.9), ...
       'MutationFcn',       util_getfield_default(optsGA,'MutationFcn',{@mutationgaussian,0.2,0.4}), ...
       'ParetoFraction',    util_getfield_default(optsGA,'ParetoFraction',0.7), ...
       'StallGenLimit',     util_getfield_default(optsGA,'StallGenLimit',150), ...
       'DistanceMeasureFcn','distancecrowding', ...
       'OutputFcn',         @(options,state,flag) ga_out_best_pen(options,state,flag, scaled, params, optsEval), ...
       'UseParallel',       usePool && util_getfield_default(optsGA,'UseParallel',true), ...
       'Display','iter','PlotFcn',[], 'FunctionTolerance',1e-5);

    %% Başlangıç Popülasyonu
    % Izgaraya hizalı ilk popülasyonu oluştur (tohumlarla birlikte).
        if isempty(util_getfield_default(optsGA,'InitialPopulationMatrix',[]))
            step_vec = [0.1 NaN 0.01 0.02 0.01 0.01 0.05 1 25 1 0.5 5 NaN 0.05 0.1];
            P0 = util_initial_pop_grid(lb, ub, options.PopulationSize, step_vec);
            % Güncel tasarım uzayı ile uyumlu tohum seti
            seed = [
                1.20  4  0.45 2.60 0.62 0.78 0.95 130 600 130 11  95  6 0.85 0.25;
                1.60  5  0.55 3.00 0.65 0.82 1.05 150 600 150 12 110  8 1.00 0.50;
                2.00  6  0.65 3.50 0.70 0.88 1.15 170 600 170 13 120 10 1.20 1.00;
                2.40  7  0.75 4.00 0.78 0.92 1.25 190 600 190 14 130 12 1.40 1.50;
                2.80  8  0.85 4.50 0.85 0.96 1.35 200 600 210 15 140 14 1.60 2.00;
            ];
            seed = max(min(seed, ub), lb);
            ns = min(size(seed,1), size(P0,1));
            P0(1:ns,:) = seed(1:ns,:);
            options = optimoptions(options,'InitialPopulationMatrix', P0);
        end

    [X,F,exitflag,output] = gamultiobj(obj, numel(lb), [],[],[],[], lb, ub, [], IntCon, options);
    gaout = struct('exitflag',exitflag,'output',output);

    %% Sonuçların Paketlenmesi
    % GA tamamlandıktan sonra sonuçları dosyalara kaydet.
    tstamp = datestr(now,'yyyymmdd_HHMMSS_FFF');
    outdir = fullfile('out', ['ga_' tstamp]);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    date_str = tstamp;
    front = struct('X',X,'F',F,'options',options,'date_str',date_str);
    save(fullfile(outdir,'ga_front.mat'), '-struct', 'front', '-v7.3');

    % === Re-evaluate Pareto designs to collect metrics per-row ===
    nF = size(X,1);
    x10_max_damperli  = zeros(nF,1);  a10abs_max_damperli = zeros(nF,1);
    dP95   = zeros(nF,1);  Qcap95 = zeros(nF,1); cav_pct = zeros(nF,1);
    T_end  = zeros(nF,1);  mu_end  = zeros(nF,1);
    PF_p95 = zeros(nF,1);  Q_q50 = zeros(nF,1);  Q_q95 = zeros(nF,1);
    dP50   = zeros(nF,1);
    energy_tot_sum   = zeros(nF,1);  E_orifice_sum = zeros(nF,1);   E_struct_sum  = zeros(nF,1); E_ratio = zeros(nF,1); P_mech_sum = zeros(nF,1);
    PFA_mean   = zeros(nF,1);  IDR_mean  = zeros(nF,1);

    % ceza bileşenleri (eval ile aynı)
    pen     = zeros(nF,1);
    pen_dP  = zeros(nF,1); pen_Qcap = zeros(nF,1); pen_cav = zeros(nF,1); pen_T = zeros(nF,1); pen_mu = zeros(nF,1);
    lambda  = util_getfield_default(optsEval,'penalty_scale',5);   % daha yumuşak ceza
    pwr     = util_getfield_default(optsEval,'penalty_power',1.0);
    W       = util_getfield_default(optsEval,'penalty_weights', struct('dP',1,'Qcap',1,'cav',1.5,'T',0.5,'mu',0.5));
    rel = @(v,lim) max(0,(v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0,(lim - v)./max(lim,eps)).^pwr;

    Opost = struct('thr', optsEval.thr);

    parfor i = 1:nF
        Xi = quant_clamp_x(X(i,:));
        Pi = decode_params_from_x(params, Xi);
        Si = run_batch_windowed(scaled, Pi, Opost);
        metrics = summarize_metrics_table(Si.table, Opost, lambda, pwr, W);

        PFA_mean(i) = metrics.PFA_mean;
        IDR_mean(i) = metrics.IDR_mean;
        x10_max_damperli(i) = metrics.x10_max_damperli;
        a10abs_max_damperli(i) = metrics.a10abs_max_damperli;
        dP95(i)   = metrics.dP95;
        Qcap95(i) = metrics.Qcap95;
        cav_pct(i) = metrics.cav_pct;
        T_end(i)  = metrics.T_end;
        mu_end(i) = metrics.mu_end;
        PF_p95(i) = metrics.PF_p95;
        Q_q50(i)  = metrics.Q_q50;
        Q_q95(i)  = metrics.Q_q95;
        dP50(i)   = metrics.dP50;
        energy_tot_sum(i) = metrics.energy_tot_sum;
        E_orifice_sum(i)  = metrics.E_orifice_sum;
        E_struct_sum(i)   = metrics.E_struct_sum;
        E_ratio(i)        = metrics.E_ratio;
        P_mech_sum(i)     = metrics.P_mech_sum;
        pen_dP(i) = metrics.pen_dP;
        pen_Qcap(i) = metrics.pen_Qcap;
        pen_cav(i) = metrics.pen_cav;
        pen_T(i)   = metrics.pen_T;
        pen_mu(i)  = metrics.pen_mu;
        pen(i)     = metrics.pen;
    end

    % Satır başına dizilerden T tablosunu oluştur
    data = [X F PFA_mean IDR_mean pen pen_dP pen_Qcap pen_cav pen_T pen_mu ...
            x10_max_damperli a10abs_max_damperli dP95 Qcap95 cav_pct T_end mu_end PF_p95 ...
            Q_q50 Q_q95 dP50 energy_tot_sum E_orifice_sum E_struct_sum E_ratio P_mech_sum];
    T = array2table(data, 'VariableNames', ...
       {'d_o_mm','n_orf','PF_tau','PF_gain','Cd0','CdInf','p_exp','Lori_mm','hA_W_perK','Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref','PF_t_on', ...
        'f1','f2','PFA_mean','IDR_mean','pen','pen_dP','pen_Qcap','pen_cav','pen_T','pen_mu', ...
        'x10_max_damperli','a10abs_max_damperli','dP95','Qcap95','cav_pct','T_end','mu_end','PF_p95', ...
        'Q_q50','Q_q95','dP50','energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio','P_mech_sum'});

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

function [f, meta, details] = eval_design_fast(x, scaled, params_base, optsEval)
    details = struct();
    % ızgaralara oturt
    x = x(:)';
    % kuantize et
    x(1) = util_quantize_step(x(1),0.05);  % d_o_mm
    x(3) = util_quantize_step(x(3),0.01);  % PF_tau
    x(4) = util_quantize_step(x(4),0.02);  % PF_gain
    x(5) = util_quantize_step(x(5),0.01);  % Cd0
    x(6) = util_quantize_step(x(6),0.01);  % CdInf
    x(7) = util_quantize_step(x(7),0.05);  % p_exp
    if numel(x) >= 8,  x(8)  = util_quantize_step(x(8),1);   end % Lori_mm
    if numel(x) >= 9,  x(9)  = util_quantize_step(x(9),25);  end % hA_W_perK
    if numel(x) >=10,  x(10) = util_quantize_step(x(10),1);  end % Dp_mm
    if numel(x) >=11,  x(11) = util_quantize_step(x(11),0.5);end % d_w_mm
    if numel(x) >=12,  x(12) = util_quantize_step(x(12),5);  end % D_m_mm
    if numel(x) >=13,  x(13) = round(x(13));                end % n_turn
    if numel(x) >=14,  x(14) = util_quantize_step(x(14),0.05); end
    if numel(x) >=15,  x(15) = util_quantize_step(x(15),0.25); end % mu_ref
    x(2)  = round(max(x(2),1));        % n_orf tam sayı, >=1
    if numel(x) >=13, x(13) = round(max(x(13),1)); end

    

    persistent memo;
    if isempty(memo), memo = containers.Map(); end
    key = jsonencode([x, dataset_signature(scaled)]);

    O = struct();
    if nargin >= 4 && ~isempty(optsEval), O = optsEval; end

    if isKey(memo, key)
        meta = memo(key);
        f = meta.f;
        if nargout > 2
            P_cached = decode_params_from_x(params_base, meta.x);
            details = run_batch_windowed(scaled, P_cached, O);
        end
        return;
    end

    P = decode_params_from_x(params_base, x);

    % Güvenli değerlendirme (GA sırasında IO yok)
    S = run_batch_windowed(scaled, P, O);

    % --- HARD FILTER (erken eleme) ---
    dP95v = S.table.dP95;
    qcapv = S.table.Qcap95;
    cavv  = S.table.cav_pct;
    if any(dP95v > 1e9) || any(qcapv > 0.90) || any(cavv > 0.01)
        f = [1e6, 1e6];
        meta = struct('x',x,'f',f,'hard_kill',true, ...
                      'pen',0,'pen_dP',0,'pen_Qcap',0,'pen_cav',0,'pen_T',0,'pen_mu',0);
        memo(key) = meta;
        if nargout > 2
            details = S;
        end
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
    lambda  = util_getfield_default(optsEval,'penalty_scale',10);
    pwr     = util_getfield_default(optsEval,'penalty_power',1.0);
    rel = @(v,lim) max(0, (v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0, (lim - v)./max(lim,eps)).^pwr;
    W = util_getfield_default(optsEval,'penalty_weights', struct('dP',1,'Qcap',1,'cav',2,'T',0.5,'mu',0.5));
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
    Penalty   = lambda*pen;
    pen_parts = struct('dP',pen_dP,'Qcap',pen_Qcap,'cav',pen_cav,'T',pen_T,'mu',pen_mu);
    if ~isfinite(Penalty), Penalty = 0; end
    pf = {'dP','Qcap','cav','T','mu'};
    for ii=1:numel(pf)
        fn = pf{ii};
        if ~isfield(pen_parts,fn) || ~isfinite(pen_parts.(fn)), pen_parts.(fn) = 0; end
    end
    x10_max_damperli_local = 0; a10abs_max_damperli_local = 0;
    if isfield(S,'table') && istable(S.table)
        if ismember('x10_max_damperli', S.table.Properties.VariableNames)
            x10_max_damperli_local = max(S.table.x10_max_damperli);
        end
        if ismember('a10abs_max_damperli', S.table.Properties.VariableNames)
            a10abs_max_damperli_local = max(S.table.a10abs_max_damperli);
        end
    end
    meta = struct('x',x,'f',f,'PFA_mean',f1,'IDR_mean',f2, ...
                 'pen',Penalty,'pen_dP',pen_parts.dP,'pen_Qcap',pen_parts.Qcap, ...
                 'pen_cav',pen_parts.cav,'pen_T',pen_parts.T,'pen_mu',pen_parts.mu, ...
                 'pen_parts',pen_parts,'x10_max_damperli',x10_max_damperli_local, ...
                 'a10abs_max_damperli',a10abs_max_damperli_local);
    % === Penaltı sürücüleri ve diagnostikleri ekle (varsa) ===
        if isfield(S,'table') && istable(S.table)
            candCols = { ...
              'dP95','Qcap95','cav_pct','T_end','mu_end', ...
              'PF_p95', ...
              'Q_q50','Q_q95','dP50', ...
              'energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio','P_mech_sum', ...
              'x10_max_damperli','a10abs_max_damperli','Re_max' ...
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
    if nargout > 2
        details = S;
    end

end

function T = prepend_baseline_row(T, params, scaled, Opost, lambda, pwr, W)
    X0 = encode_params_to_x(params);
    X0 = quant_clamp_x(X0);
    [f0, ~, S0] = eval_design_fast(X0, scaled, params, Opost);
    tbl0 = table();
    if isstruct(S0) && isfield(S0,'table') && istable(S0.table)
        tbl0 = S0.table;
    end
    metrics = summarize_metrics_table(tbl0, Opost, lambda, pwr, W);

    vn = T.Properties.VariableNames;
    T0 = T(1,:);
    T0{1,:} = nan(1,width(T0));

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
    if numel(X0) >= 15
        assign('PF_t_on', X0(15));
    end

    assign('f1', f0(1));
    assign('f2', f0(2));
    assign('PFA_mean', metrics.PFA_mean);
    assign('IDR_mean', metrics.IDR_mean);

    assign('pen',     metrics.pen);
    assign('pen_dP',  metrics.pen_dP);
    assign('pen_Qcap',metrics.pen_Qcap);
    assign('pen_cav', metrics.pen_cav);
    assign('pen_T',   metrics.pen_T);
    assign('pen_mu',  metrics.pen_mu);

    assign('x10_max_damperli', metrics.x10_max_damperli);
    assign('a10abs_max_damperli', metrics.a10abs_max_damperli);
    assign('dP95', metrics.dP95);
    assign('Qcap95', metrics.Qcap95);
    assign('cav_pct', metrics.cav_pct);
    assign('T_end', metrics.T_end);
    assign('mu_end', metrics.mu_end);
    assign('PF_p95', metrics.PF_p95);
    assign('Q_q50', metrics.Q_q50);
    assign('Q_q95', metrics.Q_q95);
    assign('dP50', metrics.dP50);
    assign('E_orifice_sum', metrics.E_orifice_sum);
    assign('E_struct_sum', metrics.E_struct_sum);
    assign('energy_tot_sum', metrics.energy_tot_sum);
    assign('E_ratio', metrics.E_ratio);
    assign('P_mech_sum', metrics.P_mech_sum);

    T = [T0; T];

    function assign(name, val)
        if ~ismember(name, vn)
            return;
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

function write_pareto_results(T, outdir)
    writetable(T, fullfile(outdir,'ga_front.csv'));
end


function metrics = summarize_metrics_table(tbl, Opost, lambda, pwr, W)
    if nargin < 1 || isempty(tbl)
        tbl = table();
    end
    if nargin < 2 || isempty(Opost)
        Opost = struct();
    end
    if nargin < 3 || isempty(lambda)
        lambda = 0;
    end
    if nargin < 4 || isempty(pwr)
        pwr = 1;
    end
    if nargin < 5 || isempty(W)
        W = struct('dP',1,'Qcap',1,'cav',1,'T',1,'mu',1);
    end
    rel = @(v,lim) max(0,(v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0,(lim - v)./max(lim,eps)).^pwr;

    metrics = struct();
    v_PFA = get_numeric('PFA', 0);
    v_IDR = get_numeric('IDR', 0);
    v_x10 = get_numeric('x10_max_damperli', 0);
    v_a10 = get_numeric('a10abs_max_damperli', 0);
    v_dP95 = get_numeric('dP95', 0);
    v_Qcap = get_numeric('Qcap95', 0);
    v_cav  = get_numeric('cav_pct', 0);
    v_T_end = get_numeric('T_end', 0);
    v_mu   = get_numeric('mu_end', 1);
    v_PF_p95 = get_numeric('PF_p95', 0);
    v_Q_q50  = get_numeric('Q_q50', 0);
    v_Q_q95  = get_numeric('Q_q95', 0);
    v_dP50   = get_numeric('dP50', 0);
    v_E_orifice = get_numeric('E_orifice_sum', 0);
    v_E_struct  = get_numeric('E_struct_sum', 0);
    v_P_mech    = get_numeric('P_mech_sum', 0);

    v_PF_p95(~isfinite(v_PF_p95)) = 0;
    v_Q_q50(~isfinite(v_Q_q50))   = 0;
    v_Q_q95(~isfinite(v_Q_q95))   = 0;
    v_dP50(~isfinite(v_dP50))     = 0;
    v_E_orifice(~isfinite(v_E_orifice)) = 0;
    v_E_struct(~isfinite(v_E_struct))   = 0;
    v_P_mech(~isfinite(v_P_mech))       = 0;

    metrics.PFA_mean = mean(v_PFA(:));
    metrics.IDR_mean = mean(v_IDR(:));
    metrics.x10_max_damperli = max(v_x10(:));
    metrics.a10abs_max_damperli = max(v_a10(:));
    metrics.dP95   = max(v_dP95(:));
    metrics.Qcap95 = max(v_Qcap(:));
    metrics.cav_pct = max(v_cav(:));
    metrics.T_end  = max(v_T_end(:));
    metrics.mu_end = min(v_mu(:));
    metrics.PF_p95 = max(v_PF_p95(:));
    metrics.Q_q50  = max(v_Q_q50(:));
    metrics.Q_q95  = max(v_Q_q95(:));
    metrics.dP50   = max(v_dP50(:));
    metrics.E_orifice_sum = sum(v_E_orifice(:));
    metrics.E_struct_sum  = sum(v_E_struct(:));
    metrics.energy_tot_sum = metrics.E_orifice_sum + metrics.E_struct_sum;
    metrics.E_ratio = 0;
    if metrics.E_struct_sum > 0
        metrics.E_ratio = metrics.E_orifice_sum / max(metrics.E_struct_sum, eps);
    end
    metrics.P_mech_sum = sum(v_P_mech(:));

    thr = util_getfield_default(Opost, 'thr', struct());
    metrics.pen_dP   = mean(rel(v_dP95(:), util_getfield_default(thr,'dP95_max',inf)));
    metrics.pen_Qcap = mean(rel(v_Qcap(:), util_getfield_default(thr,'Qcap95_max',inf)));
    cav_lim = util_getfield_default(thr,'cav_pct_max', inf);
    if cav_lim <= 0
        metrics.pen_cav = mean(max(0, v_cav(:)).^pwr);
    else
        metrics.pen_cav = mean(rel(v_cav(:), cav_lim));
    end
    metrics.pen_T    = mean(rel(v_T_end(:), util_getfield_default(thr,'T_end_max',inf)));
    metrics.pen_mu   = mean(rev(v_mu(:), util_getfield_default(thr,'mu_end_min',0)));
    metrics.pen = lambda*(W.dP*metrics.pen_dP + W.Qcap*metrics.pen_Qcap + ...
                          W.cav*metrics.pen_cav + W.T*metrics.pen_T + W.mu*metrics.pen_mu);

    function arr = get_numeric(varName, defaultVal)
        if istable(tbl) && ismember(varName, tbl.Properties.VariableNames)
            v = tbl.(varName);
            if isempty(v)
                v = defaultVal;
            end
        else
            v = defaultVal;
        end
        if ~isnumeric(v) || isempty(v)
            arr = defaultVal;
        else
            arr = v;
        end
    end
end

function sig = dataset_signature(scaled)
    if nargin < 1 || isempty(scaled)
        sig = 0;
        return;
    end
    IM = arrayfun(@(s) util_getfield_default(s,'IM',0), scaled);
    PGA = arrayfun(@(s) util_getfield_default(s,'PGA',0), scaled);
    sig = sum(double(IM(:))) + sum(double(PGA(:))) + numel(scaled);
end

function x = encode_params_to_x(params)
    orf = util_getfield_default(params,'orf', struct());
    thermal = util_getfield_default(params,'thermal', struct());
    cfg = util_getfield_default(params,'cfg', struct());
    pf = util_getfield_default(cfg,'PF', struct());

    x = nan(1,15);
    x(1) = 1e3 * util_getfield_default(orf,'d_o',3.0e-3);
    x(2) = util_getfield_default(params,'n_orf',6);
    x(3) = util_getfield_default(pf,'tau',1.0);
    x(4) = util_getfield_default(pf,'gain',0.85);
    x(5) = util_getfield_default(orf,'Cd0',0.61);
    x(6) = util_getfield_default(orf,'CdInf',0.80);
    x(7) = util_getfield_default(orf,'p_exp',1.10);
    x(8) = 1e3 * util_getfield_default(params,'Lori',0.15);
    x(9) = util_getfield_default(thermal,'hA_W_perK',600);
    x(10) = 1e3 * util_getfield_default(params,'Dp',0.150);
    x(11) = 1e3 * util_getfield_default(params,'d_w',0.012);
    x(12) = 1e3 * util_getfield_default(params,'D_m',0.120);
    x(13) = util_getfield_default(params,'n_turn',8);
    x(14) = util_getfield_default(params,'mu_ref',0.9);
    x(15) = util_getfield_default(pf,'t_on',0.75);

    if ~isfinite(x(2))
        x(2) = 6;
    end
end




function xq = quant_clamp_x(x)
    % Apply the same quantization/clamps as in eval
    xq = x;
    if isvector(x), xq = x(:)'; end
    xq(:,1) = util_quantize_step(xq(:,1),0.05);
    xq(:,3) = util_quantize_step(xq(:,3),0.01);
    xq(:,4) = util_quantize_step(xq(:,4),0.02);
    xq(:,5) = util_quantize_step(xq(:,5),0.01);
    xq(:,6) = util_quantize_step(xq(:,6),0.01);
    xq(:,7) = util_quantize_step(xq(:,7),0.05);
    if size(xq,2) >= 8,  xq(:,8)  = util_quantize_step(xq(:,8),1);   end
    if size(xq,2) >= 9,  xq(:,9)  = util_quantize_step(xq(:,9),25);  end
    if size(xq,2) >=10,  xq(:,10) = util_quantize_step(xq(:,10),1);  end
    if size(xq,2) >=11,  xq(:,11) = util_quantize_step(xq(:,11),0.5);end
    if size(xq,2) >=12,  xq(:,12) = util_quantize_step(xq(:,12),5);  end
    if size(xq,2) >=13,  xq(:,13) = round(xq(:,13));                end
    if size(xq,2) >=14,  xq(:,14) = util_quantize_step(xq(:,14),0.05); end
    if size(xq,2) >=15,  xq(:,15) = util_quantize_step(xq(:,15),0.25); end
    xq(:,2) = round(max(xq(:,2),1));
    if size(xq,2) >=13, xq(:,13) = round(max(xq(:,13),1)); end
    if isvector(x), xq = xq(:)'; end
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

function [scaled, params, T1] = prepare_inputs(optsGA)
%PREPARE_INPUTS Self-contained input assembly for the GA driver.
%   [SCALED, PARAMS, T1] = PREPARE_INPUTS(OPTSGA) builds the baseline
%   parameter struct and loads the ground-motion set using the options in
%   OPTSGA. The helper uses the local PARAMETRELER_DEFAULTS and
%   LOAD_GROUND_MOTIONS functions so no external scripts are required.

if nargin < 1 || isempty(optsGA)
    optsGA = struct();
end

[params, T1] = parametreler_defaults();
load_opts = util_getfield_default(optsGA, 'load_opts', []);
if isempty(load_opts)
    [~, scaled] = load_ground_motions(T1);
else
    [~, scaled] = load_ground_motions(T1, load_opts);
end

end

function [params, T1] = parametreler_defaults()
%PARAMETRELER_DEFAULTS Baseline structural and damper parameters.
%   [PARAMS, T1] = PARAMETRELER_DEFAULTS() reproduces the configuration from
%   the original parametreler.m script and returns the assembled parameter
%   struct together with the first-mode period T1.

%% --- 1) Structure (10 story shear building) ---
n  = 10;
m  = 360e3 * ones(n,1);
k  = 6.5e8 * ones(n,1);
c  = 6.2e6 * ones(n,1);
story_height = 3.0;

M  = diag(m);
K  = zeros(n);
C0 = zeros(n);
for i = 1:n
    kL = k(i);    cL = c(i);
    if i < n
        kU = k(i+1); cU = c(i+1);
    else
        kU = 0;    cU = 0;
    end
    K(i,i)   = kL + kU;   C0(i,i)   = cL + cU;
    if i > 1
        K(i,i-1) = -kL;  C0(i,i-1) = -cL;
    end
    if i < n
        K(i,i+1) = -kU;  C0(i,i+1) = -cU;
    end
end

[~,D] = eig(K,M);
w  = sqrt(sort(diag(D),'ascend'));
T1 = 2*pi/w(1);

%% --- 2) Damper geometry and material data ---
Dp     = 0.125;
Lgap   = 0.055;
d_o    = 3.0e-3;
Lori   = 0.10;
mu_ref = 0.9;

Kd     = 1.6e9;
Ebody  = 2.1e11;
Gsh    = 79e9;
d_w    = 12e-3;
D_m    = 80e-3;
n_turn = 8;

%% --- 3) Orifice and thermal parameters ---
rho   = 850;
n_orf = 6;

orf = struct('Cd0', 0.61, ...
             'CdInf', 0.80, ...
             'Rec', 3000, ...
             'p_exp', 1.1, ...
             'p_amb', 1.0e5, ...
             'p_cav_eff', 2.0e3, ...
             'cav_sf', 0.90, ...
             'd_o', d_o, ...
             'veps', 0.10);

T0_C    = 25;
T_ref_C = 25;
b_mu    = -0.013;
thermal = struct();
thermal.hA_W_perK = 450;
thermal.T_env_C   = 25;
thermal.max_iter  = 3;
thermal.tol_K     = 0.5;
thermal.relax     = 0.5;
thermal.dT_max    = 80;

steel_to_oil_mass_ratio = 1.5;
n_dampers_per_story    = 1;
story_mask             = ones(n-1,1);
cp_oil   = 1800;
cp_steel = 500;
resFactor = 12;

c_lam_cap      = 2e7;
c_lam_min_frac = 0.05;
c_lam_min_abs  = 1e5;

cfg = struct();
cfg.PF = struct('mode','ramp','tau',1.0,'gain',0.85,'t_on',0,'auto_t_on',true,'k',0.01);
cfg.on = struct('pressure_force', true, 'mu_floor', false, 'pf_resistive_only', false);
cfg.compat_simple = false;
cfg.num = struct('softmin_eps', 1e4, 'mu_min_phys', 0.6, 'dP_cap', NaN);

params = struct('M',M,'C0',C0,'K',K,'Dp',Dp,'Lgap',Lgap,'d_w',d_w, ...
    'D_m',D_m,'n_turn',n_turn,'mu_ref',mu_ref,'Lori',Lori,'Kd',Kd, ...
    'Ebody',Ebody,'Gsh',Gsh,'rho',rho,'n_orf',n_orf,'orf',orf, ...
    'thermal',thermal,'T0_C',T0_C,'T_ref_C',T_ref_C,'b_mu',b_mu, ...
    'c_lam_cap',c_lam_cap,'c_lam_min_frac',c_lam_min_frac, ...
    'c_lam_min_abs',c_lam_min_abs,'cp_oil',cp_oil,'cp_steel',cp_steel, ...
    'steel_to_oil_mass_ratio',steel_to_oil_mass_ratio, ...
    'n_dampers_per_story',n_dampers_per_story,'story_mask',story_mask, ...
    'resFactor',resFactor,'cfg',cfg,'story_height',story_height);

params = build_params(params);

end

function [summary, all_out] = run_batch_windowed(scaled, params, opts)
%RUN_BATCH_WINDOWED Pencereli metriklerle birden fazla kaydı analiz eder.
%   [SUMMARY, ALL_OUT] = RUN_BATCH_WINDOWED(SCALED, PARAMS, OPTS) fonksiyonu,
%   SCALED yapı dizisindeki her yer hareketi kaydını
%   RUN_ONE_RECORD_WINDOWED ile işler ve temel metriklerin özet tablosunu
%   döndürür. ALL_OUT hücre dizisi her kayıt için tam çıktıları içerir.
%   PARAMS, yapısal ve damper özelliklerini; OPTS ise
%   RUN_ONE_RECORD_WINDOWED'e iletilen ayarları barındırır.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'thr'), opts.thr = struct(); end
opts.thr = util_default_qc_thresholds(opts.thr);

% İstenirse kayıtların yürütülme sırası yeniden düzenlenir
if isfield(opts,'order_perm')
    scaled = scaled(opts.order_perm);
end

assert(isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), ...
    'run_batch_windowed: params.thermal.hA_W_perK eksik');

%% Girdi Hazırlığı
n = numel(scaled);
vars = rbw_prepare_inputs(n, params, opts);

%% Kayıt Döngüsü
vars = rbw_record_loop(scaled, params, opts, vars);

%% Özet Tablo
summary = rbw_build_summary_table(vars, opts);
all_out = vars.all_out;

end

function vars = rbw_prepare_inputs(n, params, opts)
%PREPARE_INPUTS Çalışma için gerekli dizileri hazırla
vars = struct();
vars.all_out = cell(n,1);
vars.names    = cell(n,1);
vars.scale    = zeros(n,1);
vars.SaT1     = zeros(n,1);
vars.t5       = zeros(n,1);
vars.t95      = zeros(n,1);
vars.coverage = zeros(n,1);

policy_val = util_getfield_default(opts,'thermal_reset','each');
order_val  = util_getfield_default(opts,'order','natural');
vars.policy_col = repmat({policy_val}, n,1);
vars.order_col  = repmat({order_val}, n,1);
if isfield(opts,'cooldown_s')
    cooldown_val = opts.cooldown_s;
else
    cooldown_val = NaN;
end
vars.cooldown_col = repmat(cooldown_val, n,1);

vars.PFA    = zeros(n,1);
vars.IDR    = zeros(n,1);
vars.dP95   = zeros(n,1);
vars.Qcap95 = zeros(n,1);
vars.cav_pct = zeros(n,1);
vars.PF_p95 = zeros(n,1);
vars.zeta1_hot       = zeros(n,1);
vars.z2_over_z1_hot  = zeros(n,1);
vars.P_mech_sum      = zeros(n,1);
vars.Re_max      = zeros(n,1);
vars.Q_q95  = zeros(n,1);
vars.Q_q50  = zeros(n,1);
vars.dP50   = zeros(n,1);
vars.x10_max_damperli = zeros(n,1);
vars.a10abs_max_damperli = zeros(n,1);
vars.E_orifice_sum = zeros(n,1);
vars.E_struct_sum  = zeros(n,1);
vars.energy_tot_sum = zeros(n,1);
vars.E_ratio   = zeros(n,1);
vars.qc_pass   = false(n,1);

vars.T_start    = zeros(n,1);
vars.T_end      = zeros(n,1);
vars.mu_end     = zeros(n,1);
vars.clamp_hits = zeros(n,1);

vars.Dp_mm_col   = repmat(util_getfield_default(params,'Dp_mm',NaN), n,1);
vars.mu_ref_col  = repmat(util_getfield_default(params,'mu_ref',NaN), n,1);

end

function vars = rbw_record_loop(scaled, params, opts, vars)
%RECORD_LOOP Her kaydı pencere analizine tabi tutar
prev_ts = [];
for k = 1:numel(scaled)
    rec = scaled(k);
    out = run_one_record_windowed(rec, params, opts, prev_ts);
    prev_ts = out.ts;
    vars.all_out{k} = out; %#ok<AGROW>

    vars.names{k}    = out.name;
    vars.scale(k)    = out.scale;
    vars.SaT1(k)     = out.SaT1;
    vars.t5(k)       = out.win.t5;
    vars.t95(k)      = out.win.t95;
    vars.coverage(k) = out.win.coverage;

    vars.T_start(k)    = out.T_start;
    vars.T_end(k)      = out.T_end;
    vars.mu_end(k)     = out.mu_end;
    vars.clamp_hits(k) = out.clamp_hits;
    m_nom = out.metr;
    vars.PFA(k)    = m_nom.PFA;
    vars.IDR(k)    = m_nom.IDR;
    vars.dP95(k)   = m_nom.dP95;
    vars.Qcap95(k) = m_nom.Qcap95;
    vars.cav_pct(k)= m_nom.cav_pct;
    vars.PF_p95(k) = util_getfield_default(m_nom,'PF_p95',NaN);
    vars.zeta1_hot(k)       = util_getfield_default(m_nom,'zeta1_hot',NaN);
    vars.z2_over_z1_hot(k)  = util_getfield_default(m_nom,'z2_over_z1_hot',NaN);
    vars.P_mech_sum(k)      = util_getfield_default(m_nom,'P_mech_sum',NaN);
    vars.Re_max(k)          = util_getfield_default(m_nom,'Re_max',NaN);
    vars.Q_q95(k)  = util_getfield_default(m_nom,'Q_q95',NaN);
    vars.Q_q50(k)  = util_getfield_default(m_nom,'Q_q50',NaN);
    vars.dP50(k)   = util_getfield_default(m_nom,'dP50',NaN);
    vars.x10_max_damperli(k) = util_getfield_default(m_nom,'x10_max_damperli',NaN);
    vars.a10abs_max_damperli(k) = util_getfield_default(m_nom,'a10abs_max_damperli',NaN);
    vars.E_orifice_sum(k) = util_getfield_default(m_nom,'E_orifice_sum',NaN);
    vars.E_struct_sum(k)  = util_getfield_default(m_nom,'E_struct_sum',NaN);
    vars.energy_tot_sum(k) = util_getfield_default(m_nom,'energy_tot_sum',NaN);
    vars.E_ratio(k)   = util_getfield_default(m_nom,'E_ratio',NaN);
    vars.qc_pass(k)   = out.qc_pass;

end
end

function summary = rbw_build_summary_table(vars, opts)
%BUILD_SUMMARY_TABLE Hesaplanan metrikleri tabloya dönüştür ve QC uygula
summary = struct();

thr = opts.thr;
ok_T    = vars.T_end   <= thr.T_end_max;
ok_mu   = vars.mu_end  >= thr.mu_end_min;
ok_dP   = vars.dP95    <= thr.dP95_max;
ok_Qcap = vars.Qcap95  <  thr.Qcap95_max;
ok_cav  = vars.cav_pct == 0;
qc_reason = strings(numel(vars.names),1);
for r = 1:numel(vars.names)
    bad = {};
    if ~ok_T(r),    bad{end+1}='T';  end %#ok<AGROW>
    if ~ok_mu(r),   bad{end+1}='mu'; end %#ok<AGROW>
    if ~ok_dP(r),   bad{end+1}='dP'; end %#ok<AGROW>
    if ~ok_Qcap(r), bad{end+1}='Qcap'; end %#ok<AGROW>
    if ~ok_cav(r),  bad{end+1}='cav'; end %#ok<AGROW>
    qc_reason(r) = strjoin(bad,',');
end

summary.table = table(vars.names, vars.scale, vars.SaT1, vars.t5, vars.t95, vars.coverage, vars.policy_col, vars.order_col, vars.cooldown_col, ...
    vars.PFA, vars.IDR, vars.dP95, vars.Qcap95, vars.cav_pct, vars.PF_p95, vars.zeta1_hot, vars.z2_over_z1_hot, vars.P_mech_sum, vars.Re_max, ...
    vars.Q_q95, vars.Q_q50, vars.dP50, vars.x10_max_damperli, vars.a10abs_max_damperli, vars.E_orifice_sum, vars.E_struct_sum, vars.energy_tot_sum, vars.E_ratio, vars.qc_pass, ...
    vars.T_start, vars.T_end, vars.mu_end, vars.clamp_hits, vars.Dp_mm_col, vars.mu_ref_col, ...
    ok_T, ok_mu, ok_dP, ok_Qcap, ok_cav, qc_reason, ...
    'VariableNames', {'name','scale','SaT1','t5','t95','coverage','policy','order','cooldown_s', ...
    'PFA','IDR','dP95','Qcap95','cav_pct','PF_p95','zeta1_hot','z2_over_z1_hot','P_mech_sum','Re_max', ...
    'Q_q95','Q_q50','dP50','x10_max_damperli','a10abs_max_damperli','E_orifice_sum','E_struct_sum','energy_tot_sum','E_ratio','qc_pass', ...
    'T_start','T_end','mu_end','clamp_hits','Dp_mm','mu_ref','ok_T','ok_mu','ok_dP','ok_Qcap','ok_cav','qc_reason'});

summary.all_out = vars.all_out;
end

%% Kayıt Analiz Fonksiyonu
function out = run_one_record_windowed(rec, params, opts, prev_ts)
%RUN_ONE_RECORD_WINDOWED Pencereli metriklerle tek yer hareketini analiz eder.
%   OUT = RUN_ONE_RECORD_WINDOWED(REC, PARAMS, OPTS, PREV_TS) ölçekli kayıt
%   REC için doğrusal olmayan sönümleyici modeli çalıştırır ve Arias yoğunluğu
%   tabanlı bir zaman penceresi içinde performans metriklerini hesaplar. PARAMS
%   yapısı, yapı ve sönümleyici özelliklerini (bkz. PARAMETRELER.M) içerir.
%   OPTS davranışı kontrol eder:
%       window        - MAKE_ARIAS_WINDOW fonksiyonuna iletilen ayarlar
%       thermal_reset - 'each', 'carry' veya 'cooldown'
%       cooldown_s    - 'cooldown' modu için soğuma süresi [s]
%   PREV_TS, kayıtlara arasında termal durumu taşımak için isteğe bağlı bir
%   zaman serisi yapısına sahiptir.
%
%   OUT yapısı; name, scale, SaT1, win, metr, ts, qc_pass, T_start,
%   T_end, mu_end, clamp_hits, PFA, IDR, dP95,
%   Qcap95, cav_pct, t5, t95 ve coverage alanlarını içerir. TS
%   yalnızca COMPUTE_METRICS_WINDOWED için gerekli zaman serilerini barındırır.
%
%   Bu fonksiyon MCK_WITH_DAMPER ve COMPUTE_METRICS_WINDOWED
%   fonksiyonlarını kullanır.

% varsayılan argümanlar
if nargin < 4, prev_ts = []; end
if nargin < 3 || isempty(opts), opts = struct(); end

if isfield(opts,'thermal_reset') && strcmpi(opts.thermal_reset,'cooldown')
    if ~isfield(opts,'cooldown_s') || isempty(opts.cooldown_s) || isnan(opts.cooldown_s)
        opts.cooldown_s = 60;
    end
    opts.cooldown_s = max(opts.cooldown_s,0);
end

% QC eşikleri (eksik alanlar varsayılanlarla doldurulur)
if ~isfield(opts,'thr'), opts.thr = struct(); end
thr = util_default_qc_thresholds(opts.thr);

assert(isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), ...
    'run_one_record_windowed: params.thermal.hA_W_perK eksik');

%% Pencere Hazırlığı
if isfield(opts,'window') && ~isempty(opts.window)
    wfields = fieldnames(opts.window);
    wargs = cell(1,2*numel(wfields));
    for k = 1:numel(wfields)
        wargs{2*k-1} = wfields{k};
        wargs{2*k}   = opts.window.(wfields{k});
    end
    win = util_make_arias_window(rec.t, rec.ag, wargs{:});
else
    win = util_make_arias_window(rec.t, rec.ag);
end

% PF auto_t_on, Arias t5'e göre belirlenir (çözücü çağrısından önce)
if isfield(params,'cfg') && isstruct(params.cfg) && ...
        isfield(params.cfg,'PF') && isstruct(params.cfg.PF)
    pf = params.cfg.PF;
    if isfield(pf,'auto_t_on') && pf.auto_t_on
        t5v = NaN;
        if isfield(win,'t5')
            t5v = win.t5;
        end
        if ~(isnumeric(t5v) && isfinite(t5v))
            idxnz = find(abs(rec.ag)>1e-6,1,'first');
            if isempty(idxnz), idxnz = 1; end
            t0 = rec.t(idxnz);
            params.cfg.PF.t_on = max(t0 + 0.5, 1.0);
        else
            params.cfg.PF.t_on = t5v + 0.5;
        end
    end
end

% ham ve ölçekli kayıt karşılaştırması kaldırıldı (üst seviyede kullanılmıyor)

%% Termal Sıfırlama
Tinit = params.T0_C;
if isfield(opts,'thermal_reset')
    mode = opts.thermal_reset;
else
    mode = 'each';
end

switch mode
    case 'each'
        Tinit = params.T0_C;
    case 'carry'
        if ~isempty(prev_ts) && isfield(prev_ts,'T_oil')
            Tinit = prev_ts.T_oil(end);
        else
            Tinit = params.T0_C;
        end
    case 'cooldown'
        if ~isempty(prev_ts) && isfield(prev_ts,'T_oil')
            Tprev = prev_ts.T_oil(end);
        else
            Tprev = params.T0_C;
        end
        if isfield(opts,'cooldown_s')
            td = opts.cooldown_s;
        else
            td = 0;
        end
        C_th = util_compute_Cth_effective(params);
        hA = params.thermal.hA_W_perK;
        Tenv = params.thermal.T_env_C;
        Tinit = Tenv + (Tprev - Tenv) * exp(-hA*td / C_th);
    otherwise
        Tinit = params.T0_C;
end

%% Damperli Çözüm
% Damperli çözüm mck_with_damper fonksiyonu üzerinden yürütülür.

mu_ref_eff   = params.mu_ref;
c_lam0_eff   = params.c_lam0;
if isfield(params,'cfg') && isstruct(params.cfg) && isfield(params.cfg,'on') && isstruct(params.cfg.on) && ...
        isfield(params.cfg.on,'mu_floor') && params.cfg.on.mu_floor
    mu_min_phys = NaN;
    if isfield(params.cfg,'num') && isstruct(params.cfg.num) && ...
            isfield(params.cfg.num,'mu_min_phys') && isfinite(params.cfg.num.mu_min_phys)
        mu_min_phys = params.cfg.num.mu_min_phys;
    end
    if isfinite(mu_min_phys) && (mu_ref_eff < mu_min_phys)
        mu_ref_eff = mu_min_phys;
        scale_mu = mu_ref_eff / max(params.mu_ref, eps);
        c_lam0_eff = params.c_lam0 * scale_mu;
    end
end
    [x,a_rel,ts] = mck_with_damper(rec.t, rec.ag, params.M, params.C0, params.K, ...
    params.k_sd, c_lam0_eff, params.Lori, params.orf, params.rho, params.Ap, ...
    params.Ao, params.Qcap_big, mu_ref_eff, params.thermal, Tinit, ...
    params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, params.Lgap, ...
    params.cp_oil, params.cp_steel, params.steel_to_oil_mass_ratio, ...
    params.story_mask, params.n_dampers_per_story, params.resFactor, params.cfg);

metr = compute_metrics_windowed(rec.t, x, a_rel, rec.ag, ts, params.story_height, win, params);

% Nominal koşu için son termal ve viskozite durumları kaydedilir
T_start = Tinit;
if isfield(ts,'T_oil')
    T_end = ts.T_oil(end);
    Tmax = params.T0_C + params.thermal.dT_max;
    clamp_hits = sum(diff(ts.T_oil >= Tmax) > 0);
else
    T_end = NaN;
    clamp_hits = NaN;
end
if isfield(ts,'mu')
    mu_end = ts.mu(end);
else
    mu_end = NaN;
end

qc_pass = (metr.cav_pct <= thr.cav_pct_max) && ...
          (metr.dP95 <= thr.dP95_max) && ...
          (metr.Qcap95 <= thr.Qcap95_max) && ...
          (T_end <= thr.T_end_max) && ...
          (mu_end >= thr.mu_end_min);

%% Çıktıların Derlenmesi
out = struct('name', rec.name, ...
    'scale', util_getfield_default(rec,'scale',1), ...
    'SaT1', util_getfield_default(rec,'IM',NaN), ...
    'win', win, 'metr', metr, 'ts', ts, 'qc_pass', qc_pass, ...
    'T_start', T_start, 'T_end', T_end, 'mu_end', mu_end, ...
    'clamp_hits', clamp_hits, ...
    'PFA', metr.PFA, 'IDR', metr.IDR, 'dP95', metr.dP95, 'Qcap95', metr.Qcap95, ...
    'cav_pct', metr.cav_pct, 't5', win.t5, 't95', win.t95, 'coverage', win.coverage);
if isfield(params,'cfg') && isstruct(params.cfg) && isfield(params.cfg,'PF') && isstruct(params.cfg.PF)
    pf = params.cfg.PF;
    if isfield(pf,'t_on'),   out.PF_t_on = pf.t_on;   end
    if isfield(pf,'tau'),    out.PF_tau = pf.tau;     end
    if isfield(pf,'gain'),   out.PF_gain = pf.gain;   end
    if isfield(pf,'mode'),   out.PF_mode = pf.mode;   end
    if isfield(pf,'auto_t_on')
        out.PF_auto_t_on = logical(pf.auto_t_on);
    end
end

% Yeni parametreleri isteğe bağlı olarak kaydet
param_fields = {'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref'};
for ii = 1:numel(param_fields)
    fn = param_fields{ii};
    if isfield(params, fn)
        out.(fn) = params.(fn);
    end
end
end

function metr = compute_metrics_windowed(t, x, a_rel, ag, ts, story_height, win, params)
%COMPUTE_METRICS_WINDOWED Zaman penceresi içindeki tepkileri hesaplar.
%   METR = COMPUTE_METRICS_WINDOWED(T,X,A_REL,AG,TS,STORY_HEIGHT,WIN,PARAMS)
%   fonksiyonu, WIN.IDX tarafından tanımlanan zaman aralığında yapının
%   performansına ilişkin metrikleri üretir. METR değişkeni tepe kat mutlak
%   ivmesi, katlar arası ötelenme oranları, orifis basınç istatistikleri,
%   enerji ölçümleri ve nihai ("sıcak") damper katsayısına bağlı modal sönüm
%   oranlarını içerir.
%
%   Girdi değişkenleri T, X, A_REL ve AG sırasıyla zaman vektörü, kat yer
%   değiştirmeleri, göreli kat ivmeleri ve yer ivmesini temsil eder. TS
%   yapısal analizin ürettiği ek zaman serilerini içerir; bu yapı sadece
%   metrik hesapları için gerekli dizi alanlarını barındırır. STORY_HEIGHT her
%   katın yüksekliğidir. WIN.IDX ilgilenilen pencereyi seçen mantıksal
%   vektördür. PARAMS yapısal ve damper parametrelerini içerir.

idx = win.idx;

%% Temel Tepki
% Bu bölüm, tepe kat ivmesi ve katlar arası göreli hareketler gibi yapısal
% yanıtın temel göstergelerini hesaplar.
% Tepe katın mutlak ivmesi
a_top_abs = a_rel(:,end) + ag(:);
metr.PFA = max(abs(a_top_abs(idx)));

% Pencere içindeki tepe kat mutlak yerdeğiştirmesi ve ivmesi
metr.x10_max_damperli    = max(abs(x(idx,end)));
metr.a10abs_max_damperli = max(abs(a_top_abs(idx)));

% Katlar arası göreli ötelenme oranının maksimumu
drift = (x(:,2:end) - x(:,1:end-1)) / story_height;
metr.IDR = max(max(abs(drift(idx,:))));

%% Kat Bazlı İstatistikler
% Her kat için yüzde 50 ve yüzde 95'lik değerler ile kavitasyon oranı
% gibi istatistikler hesaplanır.

abs_dP = abs(ts.dP_orf(idx,:));
abs_Q  = abs(ts.Q(idx,:));
Qcap_ratio = abs_Q ./ max(params.Qcap_big, eps);
abs_story_force = abs(ts.story_force(idx,:));

% Reynolds sayısının pencere içindeki maksimumu
metr.Re_max = NaN;
if all(isfield(params,{'Qcap_big','Ap','Ao','rho','orf'})) && isfield(ts,'mu') && isfield(ts,'dvel') && ...
        isstruct(params.orf) && all(isfield(params.orf,{'veps','d_o'}))
    dvel_win = ts.dvel(idx,:);
    qmag = params.Qcap_big * tanh((params.Ap/params.Qcap_big) * ...
        sqrt(dvel_win.^2 + params.orf.veps^2));
    mu_win = ts.mu(idx);
    mu_mat = repmat(mu_win,1,size(qmag,2));
    Ao = params.Ao; if numel(Ao)==1, Ao = Ao*ones(1,size(qmag,2)); end
    Ao_mat = repmat(Ao,size(qmag,1),1);
    Re = (params.rho .* qmag ./ max(Ao_mat .* mu_mat,1e-9)) .* ...
         max(params.orf.d_o,1e-9);
    metr.Re_max = max(Re,[],'all');
end

% 50. ve 95. yüzdelik değerler
dP_q50         = quantile(abs_dP, 0.50);
dP_q95         = quantile(abs_dP, 0.95);
Q_q50          = quantile(abs_Q, 0.50);
Q_q95          = quantile(abs_Q, 0.95);
Qcap95 = quantile(Qcap_ratio, 0.95);
story_force_q95= quantile(abs_story_force, 0.95);

% Her kat için ortalama kavitasyon yüzdesi
cav_mean = mean(ts.cav_mask(idx,:),1);

% Hikaye kesme kuvvetine göre kritik katın belirlenmesi
[metr.story_force_q95, metr.which_story] = max(story_force_q95);
ws = metr.which_story;

% Kritik katın istatistiklerinin saklanması
metr.dP95   = dP_q95(ws);
metr.dP50   = dP_q50(ws);
metr.Q_q95  = Q_q95(ws);
metr.Q_q50  = Q_q50(ws);
metr.Qcap95 = Qcap95(ws);
metr.cav_pct         = cav_mean(ws);

%% Enerji Hesapları
% Pencerede ve tüm süreçte biriken enerji bileşenleri değerlendirilir.
w_first = find(idx,1,'first');
w_last  = find(idx,1,'last');
i0 = max(w_first-1,1);

% Toplam süreç sonundaki orifis ve yapı enerjileri
metr.E_orifice_sum = ts.E_orf(end);
metr.E_struct_sum  = ts.E_struct(end);
metr.E_ratio   = metr.E_orifice_sum / max(metr.E_struct_sum, eps);
% Toplam enerji ve ortalama mekanik güç
metr.energy_tot_sum = metr.E_orifice_sum + metr.E_struct_sum;
metr.P_mech_sum = NaN;
if isfield(ts,'P_sum') && ~isempty(ts.P_sum)
    metr.P_mech_sum = mean(ts.P_sum(idx));
end

% Seçilen pencere içindeki enerji birikimleri
metr.E_orifice_win = ts.E_orf(w_last) - ts.E_orf(i0);
metr.E_struct_win  = ts.E_struct(w_last) - ts.E_struct(i0);
metr.E_ratio_win   = metr.E_orifice_win / max(metr.E_struct_win, eps);

%% Termal Metrikler
% Yağ viskozitesi gibi termal büyüklükler değerlendirilir.
if isfield(ts,'mu')
    metr.mu_end = ts.mu(w_last);
else
    metr.mu_end = NaN;
end

% ---------------- Sıcak viskozite ile modal sönüm -------------------
req_fields = {'M','K','C0','k_sd'};
if all(isfield(params,req_fields)) && isfield(ts,'c_lam')
    M  = params.M;  K = params.K;  C0 = params.C0;
    k_sd = params.k_sd;
    nStories = size(M,1) - 1;
    c_lam = ts.c_lam;

    Kadd = zeros(size(M));
    Cadd = zeros(size(M));
    for i=1:nStories
        idx2 = [i i+1];
        % R katsayısı lineer olarak kullanılır (R), karesi (R^2) alınmaz
        k_eq = k_sd;
        c_eq = c_lam;
        kM = k_eq * [1 -1; -1 1];
        cM = c_eq * [1 -1; -1 1];
        Kadd(idx2,idx2) = Kadd(idx2,idx2) + kM;
        Cadd(idx2,idx2) = Cadd(idx2,idx2) + cM;
    end
    Ktot = K + Kadd;
    Ctot = C0 + Cadd;
    [V,D] = eig(Ktot,M);
    [w2,ord] = sort(diag(D),'ascend');
    V = V(:,ord);
    phi1 = V(:,1); phi2 = V(:,2);
    w1 = sqrt(w2(1)); w2s = sqrt(w2(2));
    n1 = phi1.'*M*phi1; n2 = phi2.'*M*phi2;
    z1 = (phi1.'*Ctot*phi1)/(2*w1*n1);
    z2 = (phi2.'*Ctot*phi2)/(2*w2s*n2);
    metr.zeta1_hot = z1;
    metr.z2_over_z1_hot = z2 / max(z1, eps);
else
    metr.zeta1_hot = NaN;
    metr.z2_over_z1_hot = NaN;
end

% ----------------------- PF metrikleri (isteğe bağlı) -----------------
if isfield(ts,'PF') && ~isempty(ts.PF)
    PF_abs = abs(ts.PF(idx,:));
    metr.PF_p95 = quantile(PF_abs(:,ws), 0.95);
end

% Parametrelerden gerekli alanlar kayıt altına alınır
param_fields = {'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref'};
for ii = 1:numel(param_fields)
    fn = param_fields{ii};
    if isfield(params, fn)
        metr.(fn) = params.(fn);
    end
end

end

function [x,a_rel,ts] = mck_with_damper(t,ag,M,C,K, k_sd,c_lam0,Lori, orf,rho,Ap,Ao,Qcap, mu_ref, ...
    thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, story_mask, ...
    n_dampers_per_story, resFactor, cfg)
%% Girdi Parametreleri
    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0 = zeros(2*n,1);
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    % Kat vektörleri
    nStories = n-1;
    mask = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi= (mask .* ndps).';
    Nvec = 1:nStories; Mvec = 2:n;

    % Başlangıç sıcaklığı ve viskozitesi
    Tser = T0_C*ones(numel(t),1);
    mu_abs = mu_ref;
    c_lam = c_lam0;

%% ODE Çözümü
    odef = @(tt,z) [ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - dev_force(tt,z(1:n),z(n+1:end),c_lam,mu_abs) - M*r*agf(tt) ) ];
    sol  = ode15s(odef,[t(1) t(end)],z0,opts);
    z    = deval(sol,t).';
    x    = z(:,1:n); v = z(:,n+1:end);

    drift = x(:,Mvec) - x(:,Nvec);
    dvel  = v(:,Mvec) - v(:,Nvec);
    % Faz 3: Lineer parçada sadece yay (laminer PF tarafında)
    F_lin = k_sd*drift;

    % Faz 6: Qcap ölçeği ve softmin eps opsiyonu
    Qcap_eff = Qcap;
    if isfield(cfg,'num') && isfield(cfg.num,'Qcap_scale') && isfinite(cfg.num.Qcap_scale)
        Qcap_eff = max(1e-9, Qcap * cfg.num.Qcap_scale);
    end
    orf_loc = orf;
    if isfield(cfg,'num') && isfield(cfg.num,'softmin_eps') && isfinite(cfg.num.softmin_eps)
        orf_loc.softmin_eps = cfg.num.softmin_eps;
    end
    params = struct('Ap',Ap,'Qcap',Qcap_eff,'orf',orf_loc,'rho',rho,...
                    'Ao',Ao,'mu',mu_abs,'F_lin',F_lin,'Lori',Lori);
    [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force(dvel, params);
    % Ek diagnostikler: dP_kv ve dP_cav (kv ve kavitasyon limitleri)
    qmag_loc = Qcap_eff * tanh( (Ap/Qcap_eff) * sqrt(dvel.^2 + orf.veps^2) );
    Re_loc   = (rho .* qmag_loc .* max(orf.d_o,1e-9)) ./ max(Ao*mu_abs,1e-9);
    Cd_loc0  = orf.Cd0; Cd_locInf = orf.CdInf; Rec_loc = orf.Rec; pexp_loc = orf.p_exp;
    Cd_loc = Cd_locInf - (Cd_locInf - Cd_loc0) ./ (1 + (Re_loc./max(Rec_loc,1)).^pexp_loc);
    Cd_loc = max(min(Cd_loc,1.2),0.2);
    dP_kv_loc = 0.5*rho .* ( qmag_loc ./ max(Cd_loc.*Ao,1e-12) ).^2;
    p_up_loc  = orf.p_amb + abs(F_lin)./max(Ap,1e-12);
    dP_cav_loc= max( (p_up_loc - orf.p_cav_eff).*orf.cav_sf, 0 );
    F_p = F_lin + F_orf;

    dp_pf = (c_lam*dvel + (F_p - k_sd*drift)) ./ Ap;
    if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
        s = tanh(20*dvel);
        dp_pf = s .* max(0, s .* dp_pf);
    end
    w_pf_vec = util_pf_weight(t, cfg) * cfg.PF.gain;
    F_p = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;

    % Geometri ölçeklendirmesi R sadece montajda uygulanır
    F_story = F_p;
    P_visc_per = c_lam * (dvel.^2);
    P_sum = sum( (P_visc_per + P_orf_per) .* multi, 2 );
    P_orf_tot = sum(P_orf_per .* multi, 2);
    % Yapısal güç kat toplam kuvvetini kullanır; ekstra çarpan kullanılmaz
    P_struct_tot = sum(F_story .* dvel, 2);
    E_orf = cumtrapz(t, P_orf_tot);
    E_struct = cumtrapz(t, P_struct_tot);

%% Termal Hesap (Phase 7: iki-düğüm diagnostik; dinamiğe geri besleme yok)
    nDtot = sum(multi);
    V_oil_per = resFactor*(Ap*(2*Lgap));
    m_oil_tot = nDtot*(rho*V_oil_per);
    m_steel_tot = steel_to_oil_mass_ratio*m_oil_tot;
    C_oil   = max(m_oil_tot*cp_oil,   eps);
    C_steel = max(m_steel_tot*cp_steel, eps);
    T_o = Tser; T_s = T0_C*ones(numel(t),1);
    hA_os   = util_getfield_default(thermal, 'hA_os',    thermal.hA_W_perK);
    hA_o_env= util_getfield_default(thermal, 'hA_o_env', thermal.hA_W_perK);
    hA_s_env= util_getfield_default(thermal, 'hA_s_env', thermal.hA_W_perK);
    dtv = diff(t);
    for k=1:numel(t)-1
        Pk = 0.5*(P_sum(k)+P_sum(k+1));
        dT_o = ( Pk - hA_os*(T_o(k)-T_s(k)) - hA_o_env*(T_o(k)-thermal.T_env_C) ) / C_oil;
        dT_s = ( + hA_os*(T_o(k)-T_s(k)) - hA_s_env*(T_s(k)-thermal.T_env_C) ) / C_steel;
        T_o(k+1) = T_o(k) + dtv(k)*dT_o;
        T_s(k+1) = T_s(k) + dtv(k)*dT_s;
        T_o(k+1) = min(max(T_o(k+1), T0_C), T0_C + thermal.dT_max);
        T_s(k+1) = min(max(T_s(k+1), T0_C), T0_C + thermal.dT_max);
    end
    mu = mu_ref*exp(b_mu*(T_o - T_ref_C));

%% Çıktı Hesabı
    % İvme için düğüm kuvvetleri
    F = zeros(numel(t),n);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a_rel = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    ts = struct('dvel', dvel, 'story_force', F_story, 'Q', Q, ...
        'dP_orf', dP_orf, 'PF', F_p, 'cav_mask', dP_orf < 0, 'P_sum', P_sum, ...
        'E_orf', E_orf, 'E_struct', E_struct, 'T_oil', T_o, 'mu', mu, 'c_lam', c_lam);

%% İç Fonksiyonlar
    function Fd = dev_force(tt,x_,v_,c_lam_loc,mu_abs_loc)
        drift_ = x_(Mvec) - x_(Nvec);
        dvel_  = v_(Mvec) - v_(Nvec);
        % Sütun yönelimli etkin parametreler
        % Faz 3: Lineer parçada sadece yay
        F_lin_ = k_sd*drift_;
        params = struct('Ap',Ap,'Qcap',Qcap,'orf',orf,'rho',rho,...
                        'Ao',Ao,'mu',mu_abs_loc,'F_lin',F_lin_,'Lori',Lori);
        [F_orf_, ~, ~, ~] = calc_orifice_force(dvel_, params);
        dp_pf_ = (c_lam_loc*dvel_ + F_orf_) ./ Ap;
        if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
            s = tanh(20*dvel_);
            dp_pf_ = s .* max(0, s .* dp_pf_);
        end
        w_pf = util_pf_weight(tt,cfg) * cfg.PF.gain;
        F_p_ = k_sd*drift_ + (w_pf .* dp_pf_) * Ap;
        % R ölçeklendirmesi yalnızca montajda uygulanır (R*multi)
        F_story_ = F_p_;
        Fd = zeros(n,1);
        Fd(Nvec) = Fd(Nvec) - F_story_;
        Fd(Mvec) = Fd(Mvec) + F_story_;
    end
    function [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force(dvel, params)
        % Phase 6 (no p-states): smoother Cd(Re) and kv-only orifice drop.
        % Laminar viscous loss is accounted in PF via c_lam*dvel; avoid double counting here.

        % Saturated volumetric flow magnitude (stability)
        qmag = params.Qcap * tanh( (params.Ap/params.Qcap) * sqrt(dvel.^2 + params.orf.veps^2) );

        % Reynolds and discharge coefficient (clamped)
        Re   = (params.rho .* qmag .* max(params.orf.d_o,1e-9)) ./ max(params.Ao*params.mu,1e-9);
        Cd0   = params.orf.Cd0;
        CdInf = params.orf.CdInf;
        p_exp = params.orf.p_exp;
        Rec   = params.orf.Rec;
        Cd    = CdInf - (CdInf - Cd0) ./ (1 + (Re./max(Rec,1)).^p_exp);
        Cd    = max(min(Cd, 1.2), 0.2);

        % kv-only drop
        dP_kv  = 0.5*params.rho .* ( qmag ./ max(Cd.*params.Ao,1e-12) ).^2;

        % Cavitation soft-limit via softmin
        p_up   = params.orf.p_amb + abs(params.F_lin)./max(params.Ap,1e-12);
        dP_cav = max( (p_up - params.orf.p_cav_eff).*params.orf.cav_sf, 0 );
        epsm = 1e5;
        if isfield(params,'orf') && isfield(params.orf,'softmin_eps') && isfinite(params.orf.softmin_eps)
            epsm = params.orf.softmin_eps;
        end
        dP_orf = util_softmin(dP_kv, dP_cav, epsm);

        % Force sign from velocity (no p-states)
        sgn = dvel ./ sqrt(dvel.^2 + params.orf.veps^2);
        F_orf = dP_orf .* params.Ap .* sgn;

        % Diagnostics (positive)
        Q = qmag;
        P_orf_per = dP_kv .* qmag;   % avoid counting laminar twice
    end
end

function [records, scaled, meta] = load_ground_motions(T1, opts)
%LOAD_GROUND_MOTIONS Zemin hareketi kayitlarini yukler, on isler ve olcekler.
%   [RAW, SCALED, META] = LOAD_GROUND_MOTIONS(T1, OPTS) fonksiyonu
%   acc_matrix.mat dosyasindaki acc_matrix1, acc_matrix2, ... degiskenlerini
%   okuyarak her kayit icin zaman vektoru, ivme ve temel buyuklukleri
%   iceren bir yapi dizisi dondurur.
%
%   T1 girilirse kayitlarin 5%% sonumlu T1 periyodundaki veya tanimli
%   bir periyot bandindaki yapay spektral ivmeleri (PSA) hesaplanir ve
%   tum kayitlar ortak bir hedef IM degerine olceklendirilir.
%
%   OPTS alanlari:
%       hp_cut   - yuksek gecis filtresi frekansi [Hz]
%       IM_mode  - 'T1' veya 'band'
%       band_fac - T1 cevresindeki band carpani [alt ust]
%       band_N   - bant icindeki periyot sayisi
%       s_bounds - olcek katsayisi sinirlari [min max]
%       verbose  - islemler hakkinda bilgi mesajlari (true/false)
%
%   Ciktilar:
%       records - ham kayitlar
%       scaled  - olceklenmis kayitlar
%       meta    - islemlere dair ozet bilgiler

%% Girdi Parametreleri
if nargin < 2, opts = struct(); end
if ~isfield(opts,'verbose'), opts.verbose = true; end
hp_cut   = util_getfield_default(opts,'hp_cut',0.05);   % yuksek gecis [Hz]
IM_mode  = util_getfield_default(opts,'IM_mode','band');
band_fac = util_getfield_default(opts,'band_fac',[0.8 1.2]);
band_N   = util_getfield_default(opts,'band_N',21);
s_bounds = util_getfield_default(opts,'s_bounds',[0.2 2.2]);

%% MAT dosyasını Yükle
raw = load('acc_matrix.mat');
fn  = fieldnames(raw);
records = struct('name',{},'t',{},'ag',{},'dt',{},'duration',{},'PGA',{},'PGV',{},'IM',{},'scale',{});

%% Ön İşleme
for k = 1:numel(fn)
    A = raw.(fn{k});
    t  = A(:,1);   ag = A(:,2);

    % 3A) Zaman vektorunu duzenle
    [t,iu] = unique(t,'stable'); ag = ag(iu);
    dt = median(diff(t));
    assert(max(abs(diff(t) - dt)) < 1e-6, 'Zaman örnekleme aralığı düzensiz');
    t  = (t(1):dt:t(end)).';
    ag = interp1(A(:,1), A(:,2), t, 'linear');

    % 3B) Temel birim ve filtreleme kontrolleri
    assert(max(abs(ag)) < 100, 'İvme büyüklüğü birim hatasına işaret ediyor');
    ag = detrend(ag,0);             % ortalama
    ag = detrend(ag,1);             % dogrusal trend
    try
        Fs = 1/dt;
        Wn = hp_cut/(Fs/2);
        [b,a] = butter(2,Wn,'high');
        ag = filtfilt(b,a,ag);
    catch
        % Butter/Filtfilt yoksa yalnizca detrend uygulanir
    end

    % 3C) Metaveri
    duration = t(end) - t(1);
    v  = cumtrapz(t,ag);
    PGA = max(abs(ag));
    PGV = max(abs(v));

    records(end+1) = struct('name',fn{k},'t',t,'ag',ag,'dt',dt, ...
                             'duration',duration,'PGA',PGA,'PGV',PGV, ...
                             'IM',[],'scale',1); %#ok<AGROW>
end

%% Yükleme Özeti
if opts.verbose, fprintf('Toplam %d zemin hareketi kaydı yüklendi:\n', numel(records)); end
for k = 1:numel(records)
    r = records(k);
    if opts.verbose, fprintf('%2d) %-12s dt=%6.4f s dur=%6.2f s PGA=%7.3f PGV=%7.3f\n', ...
        k, r.name, r.dt, r.duration, r.PGA, r.PGV); end
end

scaled = [];
meta = struct();

%% IM Hesabı ve Ölçekleme (T1 sağlandığında)
if nargin >= 1 && ~isempty(T1)
    %% IM Hesabı
    for k = 1:numel(records)
        records(k).IM = compute_IM(records(k).t, records(k).ag, IM_mode, T1, band_fac, band_N);
    end

    %% Hedef IM Seçimi ve TRIM
    IM = [records.IM];
    IM_low  = max(s_bounds(1)*IM);
    IM_high = min(s_bounds(2)*IM);
    targetIM0 = median(IM);
    targetIM  = min(max(targetIM0, IM_low), IM_high);
    doClip = (IM_low > IM_high);

    max_iter = 3; dropped = {};
    for it = 1:max_iter
        IM_low  = max(s_bounds(1)*IM);
        IM_high = min(s_bounds(2)*IM);
        if IM_low <= IM_high, break; end
        [~,idx] = max(abs(log(IM) - median(log(IM))));
        dropped{end+1} = records(idx).name; %#ok<AGROW>
        records(idx) = [];
        IM(idx) = [];
    end
    if ~isempty(dropped) && opts.verbose, fprintf('TRIM: ayıklanan uç değerler = %s\n', strjoin(dropped,', ')); end
    IM_low  = max(s_bounds(1)*IM);
    IM_high = min(s_bounds(2)*IM);
    targetIM0 = median(IM);
    targetIM  = min(max(targetIM0, IM_low), IM_high);
    doClip = (IM_low > IM_high);

    %% Ölçekleme
    scaled = records; n_clipped = 0; s_all = zeros(1,numel(records));
    for k = 1:numel(records)
        s_raw = targetIM / records(k).IM;
        if doClip
            s = min(max(s_raw, s_bounds(1)), s_bounds(2));
            if abs(s - s_raw) > 1e-12, n_clipped = n_clipped + 1; end
        else
            s = s_raw;
        end
        s_all(k) = s;

        scaled(k).ag    = s * records(k).ag;
        scaled(k).PGA   = s * records(k).PGA;
        scaled(k).PGV   = s * records(k).PGV;
        scaled(k).scale = s;
        scaled(k).s_clipped = doClip && (abs(s - s_raw) > 1e-12);
        scaled(k).trimmed   = false;

        scaled(k).IM = compute_IM(scaled(k).t, scaled(k).ag, IM_mode, T1, band_fac, band_N);
    end

    %% Hata ve Log
    err = abs([scaled.IM] - targetIM) / max(targetIM, eps) * 100;
    if strcmpi(IM_mode,'band')
        modeStr = 'band';
    else
        modeStr = 'PSA@T1';
    end
    clipCount = n_clipped * doClip;
    if opts.verbose, fprintf('Hedef IM = %.3f (%s). Maks hata = %.2f%% | uygun aralık=[%.3f, %.3f] | s_min=%.2f s_max=%.2f | KIRPILAN=%d\n', ...
        targetIM, modeStr, max(err), IM_low, IM_high, min(s_all), max(s_all), clipCount); end

    meta = struct('IM_mode', IM_mode, 'band_fac', band_fac, 's_bounds', s_bounds);
    if exist('dropped','var'), meta.TRIM_names = dropped; else, meta.TRIM_names = {}; end
end
end

%% ==== Yerel Fonksiyonlar ====
function IM = compute_IM(t, ag, mode, T1, band_fac, band_N)
%COMPUTE_IM Kayit icin hedef IM degerini hesaplar.
%   IM = COMPUTE_IM(T, AG, MODE, T1, BAND_FAC, BAND_N) fonksiyonu,
%   verilen zaman vektoru t ve ivme ag icin, secilen modda
%   yapay spektral ivme (IM) dondurur.

zeta = 0.05;
if strcmpi(mode,'band')
    Tgrid = linspace(band_fac(1)*T1, band_fac(2)*T1, band_N);
    Sa = zeros(size(Tgrid));
    for i = 1:numel(Tgrid)
        Sa(i) = calc_psa(t, ag, Tgrid(i), zeta);
    end
    IM = exp(mean(log(Sa + eps)));
else
    IM = calc_psa(t, ag, T1, zeta);
end
end

function Sa = calc_psa(t, ag, T, zeta)
%CALC_PSA Tek bir kayit icin yapay spektral ivme hesabi.
%   Sa = CALC_PSA(t, ag, T, zeta) fonksiyonu, T periyotlu ve zeta
%   sonumlu tek serbestlik dereceli osilatorun mutlak ivme tepkisinin
%   maksimum degerini dondurur.

w = 2*pi / T;
agf = griddedInterpolant(t, ag, 'linear', 'nearest');

odef = @(tt, y)[ y(2);
                 -2*zeta*w*y(2) - w*w*y(1) - agf(tt) ];

y0 = [0;0];
[~, y] = ode45(odef, t, y0);

x  = y(:,1);
xd = y(:,2);
xdd = -2*zeta*w*xd - w*w*x - ag;
abs_acc = xdd + ag;
Sa = max(abs(abs_acc));
end

function params = build_params(params)
%BUILD_PARAMS Compute derived damper and hydraulic fields once.
%   PARAMS = BUILD_PARAMS(PARAMS) fills in fields such as Ap, k_sd,
%   c_lam0, Qcap_big and c_lam_min based on the fundamental geometry and
%   material properties stored in PARAMS. The input struct is returned with
%   the additional fields populated.

if nargin < 1 || ~isstruct(params)
    params = struct();
end

% Reuse existing utility for core damper quantities
params = util_recompute_damper_params(params);

% Large orifice flow cap (per damper, adjusted for parallels)
    if isfield(params,'orf') && isfield(params.orf,'CdInf') && ...
            isfield(params,'Ao') && isfield(params,'rho')
    params.Qcap_big = max(params.orf.CdInf * params.Ao, 1e-9) * ...
        sqrt(2*1.0e9 / params.rho);
end

% Minimum laminar damping based on c_lam0
if isfield(params,'c_lam_min_abs') && isfield(params,'c_lam_min_frac') && ...
        isfield(params,'c_lam0')
    params.c_lam_min = max(params.c_lam_min_abs, ...
        params.c_lam_min_frac * params.c_lam0);
end

end

function status = parpool_hard_reset(nWorkers)
% PARPOOL_HARD_RESET parpool'u sıfırlayıp iş parçacıklarını kısıtlar.
% Eski işleri temizleyerek güvenli bir havuz açılışı sağlar.

status = true;

if nargin<1 || isempty(nWorkers)
    nWorkers = feature('numcores');
else
    validateattributes(nWorkers, {'numeric'}, {'scalar','integer','positive'});
end

try
    c = parcluster('Processes');

    if ~isempty(c.Jobs)
        delete(c.Jobs);  % çökmüş işleri temizle
    end

    p = gcp('nocreate');
    if isempty(p) || ~isvalid(p)
        parpool(c, min(nWorkers, c.NumWorkers));
    end
catch ME
    warning('parpool_hard_reset:Failed', 'Parallel pool could not be opened: %s', ME.message);
    status = false;
    return;
end

pctRunOnAll maxNumCompThreads(1);          % CPU aşırı kullanımını engelle
pctRunOnAll set(0,'DefaultFigureVisible','off');
pctRunOnAll setenv('OMP_NUM_THREADS','1');
pctRunOnAll setenv('MKL_NUM_THREADS','1');

end

function y = util_softmin(a, b, epsm)
%UTIL_SOFTMIN Smooth minimum used for cavitation blending.
y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
end

function w = util_pf_weight(t, cfg)
%UTIL_PF_WEIGHT Compute the ramp weight for pressure-force contribution.
if nargin < 2 || ~isstruct(cfg), cfg = struct(); end
if ~isfield(cfg,'on') || ~isstruct(cfg.on), cfg.on = struct(); end
if ~isfield(cfg.on,'pressure_force'), cfg.on.pressure_force = true; end
if ~isfield(cfg,'PF') || ~isstruct(cfg.PF), cfg.PF = struct(); end
if ~isfield(cfg.PF,'t_on'), cfg.PF.t_on = 0; end
if ~isfield(cfg.PF,'tau'),  cfg.PF.tau  = 1.0; end
if ~isfield(cfg,'compat_simple'), cfg.compat_simple = true; end

t = double(t);
tau_floor = 1e-6;
if cfg.compat_simple
    dt  = max(t - cfg.PF.t_on, 0);
    tau = max(cfg.PF.tau, tau_floor);
    w_local = 1 - exp(-dt ./ tau);
else
    k = util_getfield_default(cfg.PF,'k',0.01);
    k = max(k, tau_floor);
    sp_dt = (log1p(exp(-abs((t - cfg.PF.t_on)./k))) + max((t - cfg.PF.t_on)./k, 0));
    dt  = sp_dt .* k;
    sp_tau = (log1p(exp(-abs((cfg.PF.tau - tau_floor)./k))) + max((cfg.PF.tau - tau_floor)./k, 0));
    tau = sp_tau .* k + tau_floor;
    w_local = 1 - exp(-dt ./ tau);
end
w = cfg.on.pressure_force .* w_local;
end

function params = util_recompute_damper_params(params)
%UTIL_RECOMPUTE_DAMPER_PARAMS Refresh derived damper quantities.
if ~isstruct(params), return; end

if isfield(params,'Dp_mm'),    params.Dp    = params.Dp_mm/1000; end
if isfield(params,'d_w_mm'),   params.d_w   = params.d_w_mm/1000; end
if isfield(params,'D_m_mm'),   params.D_m   = params.D_m_mm/1000; end
if isfield(params,'Lori_mm'),  params.Lori  = params.Lori_mm/1000; end
if isfield(params,'orf') && isfield(params.orf,'d_o_mm')
    params.orf.d_o = params.orf.d_o_mm/1000;
end

req = {'Dp','d_w','D_m','n_turn','mu_ref','Lori','Lgap','Kd','Ebody','Gsh'};
if ~all(isfield(params,req)) || ~isfield(params,'orf') || ~isfield(params.orf,'d_o')
    return;
end

if ~isfield(params,'n_orf') && isfield(params,'orf') && isfield(params.orf,'n_orf')
    params.n_orf = params.orf.n_orf;
end
if ~isfield(params,'n_orf'), params.n_orf = 1; end
nd = 1;
if isfield(params,'nd') && isfinite(params.nd), nd = params.nd; end
if isfield(params,'n_parallel') && isfinite(params.n_parallel), nd = params.n_parallel; end
if isfield(params,'n_dampers_per_story')
    nds = params.n_dampers_per_story;
    if isnumeric(nds)
        if isscalar(nds), nd = nds; else, nd = max(1, round(max(nds(:)))); end
    end
end
nd = max(1, round(nd));

Ap = pi * params.Dp^2 / 4;
Ao_single = pi * params.orf.d_o^2 / 4;
Ao = params.n_orf * Ao_single;
Ap_eff = nd * Ap;
Ao_eff = nd * Ao;

rho_loc = util_getfield_default(params,'rho',850);
Lh = rho_loc * params.Lori / max(Ao^2, 1e-18);

k_h = params.Kd * Ap^2 / params.Lgap;
k_s = params.Ebody * Ap / params.Lgap;
k_hyd = 1 / (1/k_h + 1/k_s);
k_p = params.Gsh * params.d_w^4 / (8 * params.n_turn * params.D_m^3);
k_sd_simple = k_hyd + k_p;
k_sd_adv    = nd * (k_hyd + k_p);

c_lam0 = 12 * params.mu_ref * params.Lori * Ap^2 / (params.orf.d_o^4);

params.Ap = Ap;
params.Ao = Ao;
params.Ap_eff = Ap_eff;
params.Ao_eff = Ao_eff;
params.Lh = Lh;
params.k_p = k_p;
params.k_sd_simple = k_sd_simple;
params.k_sd_adv = k_sd_adv;
params.k_sd = k_sd_adv;
params.c_lam0 = c_lam0;
end

function win = util_make_arias_window(t, ag, varargin)
%UTIL_MAKE_ARIAS_WINDOW Arias-intensity-based window generator.
p = inputParser;
p.addParameter('p1',0.05,@(x)isscalar(x) && x>=0 && x<=1);
p.addParameter('p2',0.95,@(x)isscalar(x) && x>=0 && x<=1);
p.addParameter('pad',0.5,@(x)isscalar(x) && x>=0);
p.parse(varargin{:});
p1 = p.Results.p1; p2 = p.Results.p2; pad = p.Results.pad;
IA = cumtrapz(t, ag.^2);
IA_tot = IA(end);
if ~(isfinite(IA_tot)) || IA_tot <= eps
    t_start = t(1); t_end = t(end);
    idx = true(size(t));
    win = struct('t5',t_start,'t95',t_end,'pad',0, ...
                 't_start',t_start,'t_end',t_end,'idx',idx, ...
                 'coverage',1.0,'flag_low_arias',true);
    return;
end
IA_norm = IA / IA_tot;
[t_unique, iu] = unique(IA_norm);
interp_t = t(iu);
t5  = interp1(t_unique, interp_t, p1, 'linear');
t95 = interp1(t_unique, interp_t, p2, 'linear');
dur = t95 - t5;
if dur < 5, pad = 0.25; end
t_start = max(t(1),  t5  - pad);
t_end   = min(t(end), t95 + pad);
idx = (t >= t_start) & (t <= t_end);
coverage = trapz(t(idx), ag(idx).^2) / IA_tot;
flag_low_arias = coverage < 0.90;
win = struct('t5',t5,'t95',t95,'pad',pad, ...
             't_start',t_start,'t_end',t_end,'idx',idx, ...
             'coverage',coverage,'flag_low_arias',flag_low_arias);
end

function y = util_quantize_step(x, step)
%UTIL_QUANTIZE_STEP Snap values to a uniform grid.
if nargin < 2 || isempty(step)
    y = x;
    return;
end
y = step * round(x ./ step);
end

function v = util_getfield_default(S, fname, defaultVal)
%UTIL_GETFIELD_DEFAULT Safe field accessor with defaults.
if ~isstruct(S) || ~isfield(S, fname) || isempty(S.(fname))
    v = defaultVal;
    return;
end
val = S.(fname);
if isnumeric(val)
    if isempty(val) || any(~isfinite(val(:)))
        v = defaultVal;
    else
        v = val;
    end
else
    v = val;
end
end

function thr = util_default_qc_thresholds(optsThr)
%UTIL_DEFAULT_QC_THRESHOLDS Populate QC thresholds with defaults.
if nargin < 1 || isempty(optsThr)
    optsThr = struct();
end
thr = struct();
thr.dP95_max   = util_getfield_default(optsThr,'dP95_max',50e6);
thr.Qcap95_max = util_getfield_default(optsThr,'Qcap95_max',0.5);
thr.cav_pct_max= util_getfield_default(optsThr,'cav_pct_max',0);
thr.T_end_max  = util_getfield_default(optsThr,'T_end_max',75);
thr.mu_end_min = util_getfield_default(optsThr,'mu_end_min',0.5);
extra = setdiff(fieldnames(optsThr), fieldnames(thr));
for ii = 1:numel(extra)
    thr.(extra{ii}) = optsThr.(extra{ii});
end
end

function P = util_initial_pop_grid(lb, ub, N, steps)
%UTIL_INITIAL_POP_GRID Generate a lattice-aligned initial GA population.
d = numel(lb);
P = zeros(N, d);
for i = 1:d
    if ~isnan(steps(i))
        gs = steps(i);
        if gs <= 0
            gs = max(ub(i)-lb(i), eps)/10;
        end
        grid = lb(i):gs:ub(i);
        if numel(grid) < 2
            k = max(2, ceil((ub(i)-lb(i))/max(gs, eps)));
            grid = linspace(lb(i), ub(i), k);
        end
        idx = randi(numel(grid), [N,1]);
        vals = grid(idx);
    else
        vals = lb(i) + rand(N,1).*(ub(i)-lb(i));
    end
    P(:,i) = vals;
end
end

function Cth = util_compute_Cth_effective(params)
%UTIL_COMPUTE_CTH_EFFECTIVE Effective combined thermal capacity.
nStories = size(params.M,1) - 1;
mask = params.story_mask(:);
if numel(mask)==1, mask = mask*ones(nStories,1); end
ndps = params.n_dampers_per_story(:);
if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
multi = (mask .* ndps);
V_oil_per = params.resFactor * (params.Ap * (2*params.Lgap));
m_oil_tot = sum(multi) * (params.rho * V_oil_per);
m_steel_tot = params.steel_to_oil_mass_ratio * m_oil_tot;
Cth = max(m_oil_tot*params.cp_oil + m_steel_tot*params.cp_steel, eps);
end
