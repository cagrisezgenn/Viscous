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

    [X,F] = gamultiobj(@(x) obj(x), nvars, [],[],[],[], lb, ub, [], IntCon, options);

    % save a light artifact only (no heavy exports)
    ts = datestr(now,'yyyymmdd_HHMMSS'); outdir = fullfile('out',['ga_' ts]);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    T = array2table([X F], 'VariableNames', ...
        {'d_o_mm','n_orf','g_lo','g_mid','g_hi','PF_tau','PF_gain','PFA_w','IDR_w'});
    writetable(T, fullfile(outdir,'ga_front.csv'));
    save(fullfile(outdir,'ga_front.mat'),'X','F');
end

