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
%       name, scale, SaT1, win, metr, metr0, diag, flags, which_peak
%
%   This function requires MCK_WITH_DAMPER_TS and COMPUTE_METRICS_WINDOWED
%   on the MATLAB/Octave path.

% default arguments
if nargin < 5, prev_diag = []; end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'mu_factors'), opts.mu_factors = [0.75 1.00 1.25]; end
if ~isfield(opts,'mu_weights'), opts.mu_weights = [0.2 0.6 0.2]; end

if isfield(opts,'thermal_reset') && strcmpi(opts.thermal_reset,'cooldown')
    if ~isfield(opts,'cooldown_s') || isempty(opts.cooldown_s) || isnan(opts.cooldown_s)
        opts.cooldown_s = 60;
    end
    opts.cooldown_s = max(opts.cooldown_s,0);
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
    win = make_arias_window(rec.t, rec.ag, wargs{:});
else
    win = make_arias_window(rec.t, rec.ag);
end

flags = struct();
which_peak = NaN;

%% -------------------- Optional raw vs. scaled check ------------------
if isfield(opts,'report_raw_vs_scaled') && opts.report_raw_vs_scaled && ~isempty(rec_raw)
    [~,a_lin_scaled] = lin_MCK(rec.t, rec.ag, params.M, params.C0, params.K);
    [~,a_lin_raw]    = lin_MCK(rec_raw.t, rec_raw.ag, params.M, params.C0, params.K);
    a_top_scaled = a_lin_scaled(:,end) + rec.ag;
    a_top_raw    = a_lin_raw(:,end)   + rec_raw.ag;
    peak_scaled = max(abs(a_top_scaled));
    peak_raw    = max(abs(a_top_raw));
    flags.peak_scaled = peak_scaled;
    flags.peak_raw    = peak_raw;
    which_peak = 1 + (peak_raw > peak_scaled); % 1=scaled, 2=raw
end

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

%% ---------------------- Damperless solution -------------------------
[x0,a_rel0] = lin_MCK(rec.t, rec.ag, params.M, params.C0, params.K);
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
metr0 = [];
if ~isfield(opts,'store_metr0') || opts.store_metr0
    metr0 = compute_metrics_windowed(rec.t, x0, a_rel0, rec.ag, ts0, params.story_height, win, params0);
end

%% ----------------- Damper model with time series ---------------------
nMu = numel(mu_factors);
mu_results = struct('mu_factor',cell(1,nMu));

for i = 1:nMu
    f = mu_factors(i);
    mu_ref_eff   = params.mu_ref  * f;
    c_lam0_eff   = params.c_lam0  * f;
    [x,a_rel,ts,diag] = mck_with_damper_ts(rec.t, rec.ag, params.M, params.C0, params.K, ...
        params.k_sd, c_lam0_eff, opts.use_orifice, params.orf, params.rho, params.Ap, ...
        params.A_o, params.Qcap_big, mu_ref_eff, opts.use_thermal, params.thermal, Tinit, ...
        params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, params.Lgap, ...
        params.cp_oil, params.cp_steel, params.steel_to_oil_mass_ratio, params.toggle_gain, ...
        params.story_mask, params.n_dampers_per_story, params.resFactor, params.cfg);

    params_m = params; params_m.diag = diag;
    metr_i = compute_metrics_windowed(rec.t, x, a_rel, rec.ag, ts, params.story_height, win, params_m);

    qc_pass = (metr_i.cav_pct==0) && (metr_i.dP_orf_q95<=50e6) && ...
              (metr_i.Qcap_ratio_q95<0.5) && (metr_i.T_oil_end<=75) && ...
              (metr_i.mu_end>=0.5);

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
    Tmax = params.T0_C + params.thermal.dT_max;
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

% Weighted and worst-case summaries
fields = {'PFA_top','IDR_max','dP_orf_q95','Qcap_ratio_q95','T_oil_end','mu_end'};
weighted = struct();
worst = struct();
worst.which_mu = struct();
for k = 1:numel(fields)
    fn = fields{k};
    vals = arrayfun(@(s) s.metr.(fn), mu_results);
    weighted.(fn) = sum(mu_weights(:)'.*vals);
    [worst.(fn),idx] = max(vals);
    worst.which_mu.(fn) = mu_results(idx).mu_factor;
end

%% ----------------------- Assemble output ----------------------------
out = struct();
out.name  = rec.name;
out.scale = getfield(rec,'scale',1); %#ok<GFLD>
out.SaT1  = getfield(rec,'IM',NaN); %#ok<GFLD>
out.win   = win;
out.metr  = metr;
out.metr0 = metr0;
out.diag  = diag;
out.flags = flags;
out.which_peak = which_peak;
out.mu_results = mu_results;
out.weighted = weighted;
out.worst = worst;
out.qc_all_mu = all(arrayfun(@(s) s.qc.pass, mu_results));
out.T_start = T_start;
out.T_end = T_end;
out.mu_end = mu_end;
out.clamp_hits = clamp_hits;
end

%% ---------------------------------------------------------------------
function [x,a] = lin_MCK(t,ag,M,C,K)
%LIN_MCK Solve linear MCK system without dampers.
    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    odef = @(tt,z)[ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - M*r*agf(tt) ) ];
    z0 = zeros(2*n,1);
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    sol = ode15s(odef,[t(1) t(end)],z0,opts);
    z = deval(sol,t).';
    x = z(:,1:n);
    a = ( -(M\(C*z(:,n+1:end).' + K*z(:,1:n).')).' - ag.*r.' );
end
