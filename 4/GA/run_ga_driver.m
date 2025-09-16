function [X,F,gaout] = run_ga_driver(scaled, params, optsEval, optsGA)
%RUN_GA_DRIVER Hibrit GA sürücüsü
%   [X,F,GAOUT] = RUN_GA_DRIVER(SCALED, PARAMS, OPTSEVAL, OPTSGA)
% `scaled` ve `params` doğrudan argüman olarak verilebilir. Eğer boş
% bırakılırsa, gerekli veri seti ve parametreler `prepare_inputs`
% fonksiyonu ile hazırlanır. Uygunluk hesapları sırasında IO yapılmaz;
% yeniden ölçekleme beklenmez.

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
optsEval.thr = Utils.default_qc_thresholds(Utils.getfield_default(optsEval,'thr', struct()));
%% GA Kurulumu
% GA amaç fonksiyonu ve optimizasyon seçeneklerini hazırla.
    rng(42);

    % Karar vektörü: [d_o_mm, n_orf, PF_tau, PF_gain, Cd0, CdInf, p_exp, Lori_mm, hA_W_perK, Dp_mm, d_w_mm, D_m_mm, n_turn, mu_ref, PF_t_on]
    lb = [1.0, 4, 0.40, 2.5, 0.60, 0.75, 0.90, 120, 600, 120, 10,  90,  6, 0.80, 0.0];
    ub = [3.0, 8, 0.90, 5.0, 0.90, 1.00, 1.50, 200, 600, 240, 16, 160, 18, 2.00, 3.0];
    IntCon = [2 13];  % n_orf ve n_turn tam sayı

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
    lambda  = Utils.getfield_default(optsEval,'penalty_scale',5);   % daha yumuşak ceza
    pwr     = Utils.getfield_default(optsEval,'penalty_power',1.0);
    W       = Utils.getfield_default(optsEval,'penalty_weights', struct('dP',1,'Qcap',1,'cav',1.5,'T',0.5,'mu',0.5));
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

    thr = Utils.getfield_default(Opost, 'thr', struct());
    metrics.pen_dP   = mean(rel(v_dP95(:), Utils.getfield_default(thr,'dP95_max',inf)));
    metrics.pen_Qcap = mean(rel(v_Qcap(:), Utils.getfield_default(thr,'Qcap95_max',inf)));
    cav_lim = Utils.getfield_default(thr,'cav_pct_max', inf);
    if cav_lim <= 0
        metrics.pen_cav = mean(max(0, v_cav(:)).^pwr);
    else
        metrics.pen_cav = mean(rel(v_cav(:), cav_lim));
    end
    metrics.pen_T    = mean(rel(v_T_end(:), Utils.getfield_default(thr,'T_end_max',inf)));
    metrics.pen_mu   = mean(rev(v_mu(:), Utils.getfield_default(thr,'mu_end_min',0)));
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
    IM = arrayfun(@(s) Utils.getfield_default(s,'IM',0), scaled);
    PGA = arrayfun(@(s) Utils.getfield_default(s,'PGA',0), scaled);
    sig = sum(double(IM(:))) + sum(double(PGA(:))) + numel(scaled);
end

function x = encode_params_to_x(params)
    orf = Utils.getfield_default(params,'orf', struct());
    thermal = Utils.getfield_default(params,'thermal', struct());
    cfg = Utils.getfield_default(params,'cfg', struct());
    pf = Utils.getfield_default(cfg,'PF', struct());

    x = nan(1,15);
    x(1) = 1e3 * Utils.getfield_default(orf,'d_o',3.0e-3);
    x(2) = Utils.getfield_default(params,'n_orf',6);
    x(3) = Utils.getfield_default(pf,'tau',1.0);
    x(4) = Utils.getfield_default(pf,'gain',0.85);
    x(5) = Utils.getfield_default(orf,'Cd0',0.61);
    x(6) = Utils.getfield_default(orf,'CdInf',0.80);
    x(7) = Utils.getfield_default(orf,'p_exp',1.10);
    x(8) = 1e3 * Utils.getfield_default(params,'Lori',0.15);
    x(9) = Utils.getfield_default(thermal,'hA_W_perK',600);
    x(10) = 1e3 * Utils.getfield_default(params,'Dp',0.150);
    x(11) = 1e3 * Utils.getfield_default(params,'d_w',0.012);
    x(12) = 1e3 * Utils.getfield_default(params,'D_m',0.120);
    x(13) = Utils.getfield_default(params,'n_turn',8);
    x(14) = Utils.getfield_default(params,'mu_ref',0.9);
    x(15) = Utils.getfield_default(pf,'t_on',0.75);

    if ~isfinite(x(2))
        x(2) = 6;
    end
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
