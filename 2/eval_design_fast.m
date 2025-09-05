function [f, meta] = eval_design_fast(x, scaled, params0, optsEval)
%EVAL_DESIGN_FAST Objective wrapper for GA: coarse grid + memo + silent.
%   [f, meta] = eval_design_fast(x, scaled, params0, optsEval)
%   snaps x to coarse grids, runs a lightweight batch evaluation without
%   exports/plots/diary, applies mu-robust weighting, and caches results.

    % snap to grids
    x = x(:)';
    x(1) = quantize_step(x(1),0.1);   % d_o_mm
    x(3) = quantize_step(x(3),0.05);  % g_lo
    x(4) = quantize_step(x(4),0.05);  % g_mid
    x(5) = quantize_step(x(5),0.05);  % g_hi
    x(6) = quantize_step(x(6),0.1);   % PF_tau
    x(7) = quantize_step(x(7),0.02);  % PF_gain
    x(2) = round(x(2));               % n_orf integer

    % memoization
    persistent memo;
    if isempty(memo), memo = containers.Map(); end
    key = jsonencode(x);
    if isKey(memo, key)
        meta = memo(key);
        f = meta.f;
        return;
    end

    % decode params
    P = decode_params_from_x(params0, x);

    % lightweight evaluation: NO exports/plots/diary (fast)
    O = struct();
    if nargin >= 4 && ~isempty(optsEval)
        O = optsEval;
    end
    O.do_export = false;
    O.export = struct('plots', false, 'ds', inf); % not used when do_export=false
    O.quiet  = true;
    O.thermal_reset = 'each';
    O.order = 'natural';
    % physics toggles (default robust-on)
    O.use_orifice = true;
    O.use_thermal = true;
    % mu weighting if not provided
    if ~isfield(O,'mu_factors'), O.mu_factors = [0.75 1.00 1.25]; end
    if ~isfield(O,'mu_weights'), O.mu_weights = [0.2 0.6 0.2]; end

    % run
    S = run_batch_windowed(scaled, P, O);

    % objectives = weighted PFA & IDR (minimize)
    f1 = mean(S.table.PFA_w);
    f2 = mean(S.table.IDR_w);

    % QC penalties (soft)
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
end

