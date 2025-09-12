function [X,F] = run_ga_driver(scaled, params)
%RUN_GA_DRIVER Fast, coarse-grid NSGA-II driver with memoized evals.
%   [X,F] = RUN_GA_DRIVER(SCALED, PARAMS)
%   Uses gamultiobj with integer + quantized variables and disables all
%   exports/plots/diary during objective evaluations.

    rng(42);
    nvars = 7;  % [d_o_mm, n_orf, g_lo, g_mid, g_hi, PF_tau, PF_gain]
    lb = [1.5, 2, 0.5, 0.5, 0.5, 1.0, 0.40];
    ub = [4.0, 6, 4.0, 4.0, 4.0, 3.0, 0.80];
    IntCon = 2;  % n_orf

    % eval options (QC thresholds etc.)
    optsEval = struct;
    optsEval.thr = struct('dP95_max',50e6,'Qcap95_max',0.5, ...
                          'cav_pct_max',0,'T_end_max',90,'mu_end_min',0.35);

    obj = @(x) eval_design_fast(x, scaled, params, optsEval);

    options = optimoptions('gamultiobj', ...
        'UseParallel', true, ...
        'Display','iter', ...
        'PlotFcn', [], ...
        'MaxGenerations', 150, ...
        'PopulationSize', 60, ...
        'FunctionTolerance', 1e-3, ...
        'ConstraintTolerance', 1e-3);

    [X,F] = gamultiobj(obj, nvars, [],[],[],[], lb, ub, [], IntCon, options);

    % Package GA results (no simulations, no plots)
    date_str = datestr(now,'yyyy-mm-dd_HHMMSS');
    outdir = fullfile('out',['ga_' date_str]);
    if ~exist(outdir,'dir'), mkdir(outdir); end

    % quantize+clamp decisions for export consistency
    Xq = zeros(size(X));
    for i = 1:size(X,1)
        Xq(i,:) = quant_clamp_x(X(i,:));
    end

    % Save front (CSV): decisions then f1,f2
    Tfront = array2table([Xq F], 'VariableNames', ...
        {'d_o_mm','n_orf','g_lo','g_mid','g_hi','PF_tau','PF_gain','f1','f2'});
    writetable(Tfront, fullfile(outdir,'ga_front.csv'));

    % Options/meta and MAT
    opts_ga = options;
    thr = optsEval.thr;
    IM_mode = 'band'; s_bounds = [0.5 2.0]; band_fac = [0.8 1.2];
    mu_factors = [0.75 1.00 1.25]; mu_weights = [0.2 0.6 0.2];

    % Decode TOP-K representative designs (no simulation)
    K = min(10, size(Xq,1));
    if K > 0
        [~,ord] = sort(F(:,1), 'ascend');
        pick = unique(round(linspace(1, numel(ord), K)));
        idxK = ord(pick);
        params_list = cell(numel(idxK),1);
        topK_tbl = table('Size',[numel(idxK) 9], 'VariableTypes',repmat({'double'},1,9), ...
            'VariableNames',{'d_o_mm','n_orf','g_lo','g_mid','g_hi','PF_tau','PF_gain','A_o','Qcap_big'});
        for j = 1:numel(idxK)
            xj = Xq(idxK(j),:);
            Pj = decode_params_from_x(params, xj);
            params_list{j} = Pj; %#ok<AGROW>
            topK_tbl{j,1:7} = num2cell(xj);
            topK_tbl.A_o(j) = Pj.A_o;
            topK_tbl.Qcap_big(j) = Pj.Qcap_big;
        end
        writetable(topK_tbl, fullfile(outdir,'ga_topK.csv'));
    else
        params_list = {};
    end

    % MAT with meta
    save(fullfile(outdir,'ga_front.mat'), 'X','F','Xq','opts_ga','thr', ...
        'IM_mode','mu_factors','mu_weights','s_bounds','band_fac','date_str','params_list');

    % Minimal README
    readme = {
        sprintf('GA run @ %s', date_str)
        'Reproducibility meta:'
        sprintf('  IM_mode=%s, band_fac=[%.2f %.2f], s_bounds=[%.2f %.2f]', IM_mode, band_fac(1), band_fac(2), s_bounds(1), s_bounds(2))
        sprintf('  mu_factors=%s, mu_weights=%s', mat2str(mu_factors), mat2str(mu_weights))
        sprintf('  QC thr: %s', jsonencode(thr))
        sprintf('  GA opts: Pop=%d, Gen=%d, Seed=%d', opts_ga.PopulationSize, opts_ga.MaxGenerations, 42)
        'Note: No simulations are run during packaging.'
    };
    fid = fopen(fullfile(outdir,'README.txt'),'w'); if fid~=-1, fprintf(fid,'%s\n', readme{:}); fclose(fid); end
end

function [f, meta] = eval_design_fast(x, scaled, params0, optsEval)
    % snap to grids
    x = x(:)';
    % quantize
    x(1) = Utils.quantize_step(x(1),0.1);   % d_o_mm
    x(3) = Utils.quantize_step(x(3),0.05);  % g_lo
    x(4) = Utils.quantize_step(x(4),0.05);  % g_mid
    x(5) = Utils.quantize_step(x(5),0.05);  % g_hi
    x(6) = Utils.quantize_step(x(6),0.1);   % PF_tau
    x(7) = Utils.quantize_step(x(7),0.02);  % PF_gain
    x(2) = round(max(x(2),1));              % n_orf integer, >=1

    % clamp to GA bounds to be safe after quantization
    x(1) = min(max(x(1), 1.5), 4.0);
    x(2) = min(max(x(2), 2), 6);
    x(3) = min(max(x(3), 0.5), 4.0);
    x(4) = min(max(x(4), 0.5), 4.0);
    x(5) = min(max(x(5), 0.5), 4.0);
    x(6) = min(max(x(6), 1.0), 3.0);
    x(7) = min(max(x(7), 0.40), 0.80);

    persistent memo;
    if isempty(memo), memo = containers.Map(); end
    key = jsonencode(x);
    if isKey(memo, key)
        meta = memo(key); f = meta.f; return;
    end

    P = decode_params_from_x(params0, x);

    O = struct();
    if nargin >= 4 && ~isempty(optsEval), O = optsEval; end
    O.do_export = false;
    O.export = struct('plots', false, 'ds', inf);
    O.quiet  = true;
    O.thermal_reset = 'each';
    O.order = 'natural';
    O.use_orifice = true; O.use_thermal = true;
    if ~isfield(O,'mu_factors'), O.mu_factors = [0.75 1.00 1.25]; end
    if ~isfield(O,'mu_weights'), O.mu_weights = [0.2 0.6 0.2]; end

    % safe evaluation (no IO during GA)
    try
        S = run_batch_windowed(scaled, P, O);
    catch ME
        f = [1e6, 1e6];
        meta = struct('x',x,'f',f,'error',ME.message);
        return;
    end

    f1 = mean(S.table.PFA_w);
    f2 = mean(S.table.IDR_w);
    thr = O.thr;
    dP95v   = S.table.dP95_worst;    
    qcapv   = S.table.Qcap95_worst;
    cavv    = S.table.cav_pct_worst;
    Tendv   = S.table.T_end_worst;
    muendv  = S.table.mu_end_worst;
    pen = 0;
    if isfield(thr,'dP95_max')
        pen = pen + 1e3*sum(max(0, dP95v - thr.dP95_max)./max(thr.dP95_max,eps)).^2;
    end
    if isfield(thr,'Qcap95_max')
        pen = pen + 1e3*sum(max(0, qcapv - thr.Qcap95_max)./max(thr.Qcap95_max,eps)).^2;
    end
    if isfield(thr,'cav_pct_max')
        pen = pen + 1e3*sum(max(0, cavv - thr.cav_pct_max)./max(thr.cav_pct_max+eps,eps)).^2;
    end
    if isfield(thr,'T_end_max')
        pen = pen + 1e3*sum(max(0, Tendv - thr.T_end_max)./max(thr.T_end_max,eps)).^2;
    end
    if isfield(thr,'mu_end_min')
        pen = pen + 1e3*sum(max(0, thr.mu_end_min - muendv)./max(thr.mu_end_min,eps)).^2;
    end

    f = [f1+pen, f2+pen];
    meta = struct('x',x,'f',f,'PFA_w_mean',f1,'IDR_w_mean',f2,'pen',pen);
    memo(key) = meta;

    function P = decode_params_from_x(params0_, x_)
        d_o_mm = x_(1); n_orf = round(x_(2));
        g_lo = x_(3); g_mid = x_(4); g_hi = x_(5);
        PF_tau = x_(6); PF_gain = x_(7);
        P = params0_;
        P.orf.d_o = d_o_mm * 1e-3;         % mm -> m
        P.n_orf   = n_orf;
        P.A_o     = P.n_orf * (pi * P.orf.d_o^2 / 4);
        P.Qcap_big= max(P.orf.CdInf * P.A_o, 1e-9) * sqrt(2 * 1.0e9 / P.rho);
        n  = size(P.M,1); nStories = n-1;
        tg = ones(nStories,1) * g_mid;
        loN = min(3, nStories); if loN > 0, tg(1:loN) = g_lo; end
        hiN = min(2, nStories); if hiN > 0, tg(end-hiN+1:end) = g_hi; end
        P.toggle_gain = tg;
        if isfield(P,'cfg') && isfield(P.cfg,'PF')
            P.cfg.PF.tau  = PF_tau;
            P.cfg.PF.gain = PF_gain;
        end
    end
end

function xq = quant_clamp_x(x)
    % Apply the same quantization/clamps as in eval
    x = x(:)';
    x(1) = Utils.quantize_step(x(1),0.1);
    x(3) = Utils.quantize_step(x(3),0.05);
    x(4) = Utils.quantize_step(x(4),0.05);
    x(5) = Utils.quantize_step(x(5),0.05);
    x(6) = Utils.quantize_step(x(6),0.1);
    x(7) = Utils.quantize_step(x(7),0.02);
    x(2) = round(max(x(2),1));
    x(1) = min(max(x(1), 1.5), 4.0);
    x(2) = min(max(x(2), 2), 6);
    x(3) = min(max(x(3), 0.5), 4.0);
    x(4) = min(max(x(4), 0.5), 4.0);
    x(5) = min(max(x(5), 0.5), 4.0);
    x(6) = min(max(x(6), 1.0), 3.0);
    x(7) = min(max(x(7), 0.40), 0.80);
    xq = x;
end

function P = decode_params_from_x(params0_, x_)
    d_o_mm = x_(1); n_orf = round(x_(2));
    g_lo = x_(3); g_mid = x_(4); g_hi = x_(5);
    PF_tau = x_(6); PF_gain = x_(7);
    P = params0_;
    P.orf.d_o = d_o_mm * 1e-3;         % mm -> m
    P.n_orf   = n_orf;
    P.A_o     = P.n_orf * (pi * P.orf.d_o^2 / 4);
    P.Qcap_big= max(P.orf.CdInf * P.A_o, 1e-9) * sqrt(2 * 1.0e9 / P.rho);
    n  = size(P.M,1); nStories = n-1;
    tg = ones(nStories,1) * g_mid;
    loN = min(3, nStories); if loN > 0, tg(1:loN) = g_lo; end
    hiN = min(2, nStories); if hiN > 0, tg(end-hiN+1:end) = g_hi; end
    P.toggle_gain = tg;
    if isfield(P,'cfg') && isfield(P.cfg,'PF')
        P.cfg.PF.tau  = PF_tau;
        P.cfg.PF.gain = PF_gain;
    end
end
