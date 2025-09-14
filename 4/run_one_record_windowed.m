function out = run_one_record_windowed(rec, rec_raw, params, opts, prev_diag)
%RUN_ONE_RECORD_WINDOWED Analyse one ground-motion record with windowed metrics.
%   OUT = RUN_ONE_RECORD_WINDOWED(REC, REC_RAW, PARAMS, OPTS, PREV_DIAG) runs
%   the nonlinear damper model for the scaled record REC and computes
%   performance metrics restricted to an Arias intensity based time window.
%   REC_RAW is the unscaled record used for optional comparison when
%   OPTS.report_raw_vs_scaled is true. PARAMS bundles structural and damper
%   properties (see PARAMETRELER.M). OPTS controls behaviour:
%       window   - struct of options forwarded to MAKE_ARIAS_WINDOW
%       report_raw_vs_scaled - if true, a dampers-free solution of REC_RAW is
%                              also obtained for comparison
%       store_metr0 - if true, metrics for the damperless system are stored
%       thermal_reset - 'each', 'carry' or 'cooldown'
%       cooldown_s   - cooldown duration for 'cooldown' mode [s]
%   PREV_DIAG is an optional diagnostic structure from a previous call used
%   to carry thermal state between records when requested.
%
%   Output fields of OUT:
%       name, scale, SaT1, win, metr, diag, mu_results, weighted, worst, ts, ...
%       qc_all_mu, T_start, T_end, mu_end, clamp_hits, PFA_top, IDR_max, ...
%       dP_orf_q95, Qcap_ratio_q95, cav_pct, t5, t95, coverage
%
%   This function requires MCK_WITH_DAMPER_TS and COMPUTE_METRICS_WINDOWED
%   on the MATLAB/Octave path.

% default arguments
if nargin < 5, prev_diag = []; end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'mu_factors'), opts.mu_factors = 1.00; end
if ~isfield(opts,'mu_weights'), opts.mu_weights = 1; end

if isfield(opts,'thermal_reset') && strcmpi(opts.thermal_reset,'cooldown')
    if ~isfield(opts,'cooldown_s') || isempty(opts.cooldown_s) || isnan(opts.cooldown_s)
        opts.cooldown_s = 60;
    end
    opts.cooldown_s = max(opts.cooldown_s,0);
end

% QC thresholds (can be overridden via opts.thr)
thr_default = struct('dP95_max',50e6, 'Qcap95_max',0.5, 'cav_pct_max',0, ...
    'T_end_max',75, 'mu_end_min',0.5);
if isfield(opts,'thr') && ~isempty(opts.thr)
    thr_f = opts.thr;
    fns = fieldnames(thr_default);
    for ii=1:numel(fns)
        if ~isfield(thr_f,fns{ii}) || isempty(thr_f.(fns{ii}))
            thr_f.(fns{ii}) = thr_default.(fns{ii});
        end
    end
    thr = thr_f;
else
    thr = thr_default;
end

assert(numel(opts.mu_factors)==numel(opts.mu_weights), ...
    'mu_factors and mu_weights must have same length.');
mu_weights = opts.mu_weights(:);
wsum = sum(mu_weights);
assert(wsum>0,'mu_weights sum must be > 0.');
mu_weights = mu_weights/wsum;
mu_factors = opts.mu_factors(:)';

%% ----------------------- Arias intensity window ----------------------
if isfield(opts,'window') && ~isempty(opts.window)
    wfields = fieldnames(opts.window);
    wargs = cell(1,2*numel(wfields));
    for k = 1:numel(wfields)
        wargs{2*k-1} = wfields{k};
        wargs{2*k}   = opts.window.(wfields{k});
    end
    win = Utils.make_arias_window(rec.t, rec.ag, wargs{:});
else
    win = Utils.make_arias_window(rec.t, rec.ag);
end

% PF auto_t_on based on Arias t5 (before any solver call)
try
    if isfield(params,'cfg') && isfield(params.cfg,'PF') && ...
       isfield(params.cfg.PF,'auto_t_on') && params.cfg.PF.auto_t_on
        t5v = NaN; if isfield(win,'t5'), t5v = win.t5; end
        if ~(isnumeric(t5v) && isfinite(t5v))
            idxnz = find(abs(rec.ag)>1e-6,1,'first');
            if isempty(idxnz), idxnz = 1; end
            t0 = rec.t(idxnz);
            params.cfg.PF.t_on = max(t0 + 0.5, 1.0);
        else
            params.cfg.PF.t_on = t5v + 0.5;
        end
    end
catch
    % ignore errors; leave params unchanged
end

% raw-vs-scaled comparison removed (unused upstream)

%% -------------------- Thermal reset handling -------------------------
Tinit = params.T0_C;
if isfield(opts,'thermal_reset')
    mode = opts.thermal_reset;
else
    mode = 'each';
end

% compute thermal capacity for cooldown option
nStories = size(params.M,1) - 1;
Rvec = params.toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
mask = params.story_mask(:);  if numel(mask)==1, mask = mask*ones(nStories,1); end
ndps = params.n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
multi = (mask .* ndps);
V_oil_per = params.resFactor * (params.Ap * (2*params.Lgap));
 m_oil_tot = sum(multi) * (params.rho * V_oil_per);
 m_steel_tot = params.steel_to_oil_mass_ratio * m_oil_tot;
 C_th = max(m_oil_tot*params.cp_oil + m_steel_tot*params.cp_steel, eps);

switch mode
    case 'each'
        Tinit = params.T0_C;
    case 'carry'
        if ~isempty(prev_diag) && isfield(prev_diag,'T_oil')
            Tinit = prev_diag.T_oil(end);
        else
            Tinit = params.T0_C;
        end
    case 'cooldown'
        if ~isempty(prev_diag) && isfield(prev_diag,'T_oil')
            Tprev = prev_diag.T_oil(end);
        else
            Tprev = params.T0_C;
        end
        if isfield(opts,'cooldown_s')
            td = opts.cooldown_s;
        else
            td = 0;
        end
        hA = params.thermal.hA_W_perK;
        Tenv = params.thermal.T_env_C;
        Tinit = Tenv + (Tprev - Tenv) * exp(-hA*td / C_th);
    otherwise
        Tinit = params.T0_C;
end

params.thermal.T0_C = Tinit;

%% ---------------------- Damperless solution -------------------------
[x0,a_rel0] = Utils.lin_MCK(rec.t, rec.ag, params.M, params.C0, params.K);
ts0 = struct('dP_orf',zeros(numel(rec.t),nStories), ...
             'Q',zeros(numel(rec.t),nStories), ...
             'Qcap_ratio',zeros(numel(rec.t),nStories), ...
             'story_force',zeros(numel(rec.t),nStories), ...
             'cav_mask',false(numel(rec.t),nStories), ...
             'E_orf',zeros(numel(rec.t),1), ...
             'E_struct',zeros(numel(rec.t),1));
params0 = params; params0.diag = struct('T_oil',zeros(numel(rec.t),1), ...
                                       'mu',zeros(numel(rec.t),1), ...
                                       'c_lam',0);
% store_metr0 removed (unused upstream)

%% ----------------- Damper model with time series ---------------------
nMu = numel(mu_factors);
mu_results = struct('mu_factor',cell(1,nMu));

for i = 1:nMu
    f = mu_factors(i);
    mu_ref_eff   = params.mu_ref  * f;
    c_lam0_eff   = params.c_lam0  * f;
    [x,a_rel,ts,diag] = mck_with_damper_ts(rec.t, rec.ag, params.M, params.C0, params.K, ...
        params.k_sd, c_lam0_eff, params.Lori, opts.use_orifice, params.orf, params.rho, params.Ap, ...
        params.A_o, params.Qcap_big, mu_ref_eff, opts.use_thermal, params.thermal, ...
        params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, params.Lgap, ...
        params.cp_oil, params.cp_steel, params.steel_to_oil_mass_ratio, params.toggle_gain, ...
        params.story_mask, params.n_dampers_per_story, params.resFactor, params.cfg);

    params_m = params; params_m.diag = diag;
    metr_i = compute_metrics_windowed(rec.t, x, a_rel, rec.ag, ts, params.story_height, win, params_m);

    qc_pass = (metr_i.cav_pct <= thr.cav_pct_max) && ...
              (metr_i.dP_orf_q95 <= thr.dP95_max) && ...
              (metr_i.Qcap_ratio_q95 <= thr.Qcap95_max) && ...
              (metr_i.T_oil_end <= thr.T_end_max) && ...
              (metr_i.mu_end >= thr.mu_end_min);

    mu_results(i).mu_factor   = f;
    mu_results(i).mu_ref_eff  = mu_ref_eff;
    mu_results(i).c_lam0_eff  = c_lam0_eff;
    mu_results(i).metr        = metr_i;
    mu_results(i).diag        = diag;
    mu_results(i).qc.pass     = qc_pass;
end

% Nominal metrics (f=1)
[~,nom_idx] = min(abs(mu_factors-1));
metr = mu_results(nom_idx).metr;
diag = mu_results(nom_idx).diag;

% Log thermal/viscosity end states for nominal run
T_start = Tinit;
if isfield(diag,'T_oil')
    T_end = diag.T_oil(end);
    Tmax = params.thermal.T0_C + params.thermal.dT_max;
    clamp_hits = sum(diff(diag.T_oil >= Tmax) > 0);
else
    T_end = NaN;
    clamp_hits = NaN;
end
if isfield(diag,'mu')
    mu_end = diag.mu(end);
else
    mu_end = NaN;
end

% Weighted and worst-case summaries (metric-specific min/max where appropriate)
fields = {'PFA_top','IDR_max','dP_orf_q95','dP_orf_q50','Q_q95','Q_q50','Qcap_ratio_q95', ...
          'cav_pct','T_oil_end','mu_end', 'x10_max_D','a10abs_max_D', ...
          'E_orifice_full','E_struct_full','E_ratio_full','PF_p95'};
weighted = struct();
worst = struct();
worst.which_mu = struct();
for kf = 1:numel(fields)
    fn = fields{kf};
    vals = arrayfun(@(s) s.metr.(fn), mu_results);
    % weighted average for all
    weighted.(fn) = sum(mu_weights(:)'.*vals);
    % worst-case selection rule
    switch fn
        case 'mu_end'
            [worst.(fn), idx] = min(vals); % smaller viscosity is worst
        otherwise
            [worst.(fn), idx] = max(vals);
    end
    worst.which_mu.(fn) = mu_results(idx).mu_factor;
end

%% ----------------------- Assemble output ----------------------------
out = struct();
out.name  = rec.name;
out.scale = Utils.getfield_default(rec,'scale',1);
out.SaT1  = Utils.getfield_default(rec,'IM',NaN);
out.win   = win;
out.metr  = metr;
out.diag  = diag;
out.mu_results = mu_results;
out.weighted = weighted;
out.worst = worst;
out.ts = ts;
out.qc_all_mu = all(arrayfun(@(s) s.qc.pass, mu_results));
out.T_start = T_start;
out.T_end = T_end;
out.mu_end = mu_end;
out.clamp_hits = clamp_hits;
% convenience telemetry fields
out.PFA_top = metr.PFA_top;
out.IDR_max = metr.IDR_max;
out.dP_orf_q95 = metr.dP_orf_q95;
out.Qcap_ratio_q95 = metr.Qcap_ratio_q95;
out.cav_pct = metr.cav_pct;
out.t5 = win.t5; out.t95 = win.t95; out.coverage = win.coverage;
% chosen PF ramp onset (if set)
try
    if isfield(params,'cfg') && isfield(params.cfg,'PF')
        if isfield(params.cfg.PF,'t_on')
            out.PF_t_on = params.cfg.PF.t_on;
        end
        if isfield(params.cfg.PF,'tau')
            out.PF_tau = params.cfg.PF.tau;
        end
        if isfield(params.cfg.PF,'gain')
            out.PF_gain = params.cfg.PF.gain;
        end
        if isfield(params.cfg.PF,'mode')
            out.PF_mode = params.cfg.PF.mode;
        end
        if isfield(params.cfg.PF,'auto_t_on')
            out.PF_auto_t_on = logical(params.cfg.PF.auto_t_on);
        end
    end
catch
end
end


