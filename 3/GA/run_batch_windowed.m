function [summary, all_out] = run_batch_windowed(scaled, params, opts)
%RUN_BATCH_WINDOWED Analyse multiple records with windowed metrics.
%   [SUMMARY, ALL_OUT] = RUN_BATCH_WINDOWED(SCALED, PARAMS, OPTS) processes each
%   ground-motion record in the struct array SCALED using
%   RUN_ONE_RECORD_WINDOWED and returns a summary table of key metrics. The
%   cell array ALL_OUT contains the full outputs for each record. PARAMS
%   bundles structural and damper properties. OPTS are forwarded to
%   RUN_ONE_RECORD_WINDOWED.
%
%   QC logs are printed for IM consistency, low Arias coverage, physical
%   plausibility of response metrics and saturation/cavitation checks.
%   After processing all records, the worst peak floor acceleration and
%   inter-story drift ratio along with the associated \mu factor are
%   reported.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'mu_factors'), opts.mu_factors = [0.75 1.00 1.25]; end
if ~isfield(opts,'mu_weights'), opts.mu_weights = [0.2 0.6 0.2]; end

do_export = isfield(opts,'do_export') && opts.do_export;
if do_export
    if isfield(opts,'outdir')
        outdir = opts.outdir;
    else
        ts = datestr(now,'yyyymmdd_HHMMSS_FFF');
        outdir = fullfile('out', ts);
    end
    if ~exist(outdir,'dir'), mkdir(outdir); end
    if ~Utils.getfield_default(opts,'quiet',false)
        diary(fullfile(outdir,'console.log'));
    end
else
    outdir = '';
end

n = numel(scaled);

all_out = cell(n,1);

names    = cell(n,1);
scale    = zeros(n,1);
SaT1     = zeros(n,1);
t5       = zeros(n,1);
t95      = zeros(n,1);
coverage = zeros(n,1);
rank_score = nan(n,1);

% policy/order info
policy_val = Utils.getfield_default(opts,'thermal_reset','each');
order_val  = Utils.getfield_default(opts,'order','natural');
policy_col = repmat({policy_val}, n,1);
order_col  = repmat({order_val}, n,1);
if isfield(opts,'cooldown_s')
    cooldown_val = opts.cooldown_s;
else
    cooldown_val = NaN;
end
cooldown_col = repmat(cooldown_val, n,1);

PFA_nom    = zeros(n,1);
IDR_nom    = zeros(n,1);
dP95_nom   = zeros(n,1);
Qcap95_nom = zeros(n,1);
cav_nom    = zeros(n,1);

PFA_w    = zeros(n,1);
IDR_w    = zeros(n,1);
dP95_w   = zeros(n,1);
Qcap95_w = zeros(n,1);
Q_q95_w  = zeros(n,1);
Q_q50_w  = zeros(n,1);
dP50_w   = zeros(n,1);

PFA_worst    = zeros(n,1);
IDR_worst    = zeros(n,1);
dP95_worst   = zeros(n,1);
Qcap95_worst = zeros(n,1);
Q_q95_worst  = zeros(n,1);
Q_q50_worst  = zeros(n,1);
dP_orf_q50_worst = zeros(n,1);
PF_p95_worst = zeros(n,1);
which_mu_PFA = zeros(n,1);
which_mu_IDR = zeros(n,1);
T_end_worst  = zeros(n,1);
mu_end_worst = zeros(n,1);
cav_pct_worst = zeros(n,1);
% New worst-case peak metrics at top story (damperli)
x10_max_D_worst = zeros(n,1);
a10abs_max_D_worst = zeros(n,1);
% Energy summaries (worst across mu)
E_orifice_sum = zeros(n,1);
E_struct_sum  = zeros(n,1);
E_ratio       = zeros(n,1);
qc_all_mu    = false(n,1);

T_start    = zeros(n,1);
T_end      = zeros(n,1);
mu_end     = zeros(n,1);
clamp_hits = zeros(n,1);

worstPFA = -inf; worstPFA_name = ''; worstPFA_mu = NaN;
worstIDR = -inf; worstIDR_name = ''; worstIDR_mu = NaN;

prev_diag = [];
for k = 1:n
    rec = scaled(k);
    out = run_one_record_windowed(rec, params, opts, prev_diag);
    prev_diag = out.diag;
    all_out{k} = out; %#ok<AGROW>

    names{k}    = out.name;
    scale(k)    = out.scale;
    SaT1(k)     = out.SaT1;
    t5(k)       = out.win.t5;
    t95(k)      = out.win.t95;
    coverage(k) = out.win.coverage;
    % rank score computed only for order='worst_first'
    if strcmpi(order_val,'worst_first')
        if isfield(out,'metr') && isfield(out.metr,'E_orifice_win')
            rank_score(k) = out.metr.E_orifice_win;
        end
    end

    T_start(k)    = out.T_start;
    T_end(k)      = out.T_end;
    mu_end(k)     = out.mu_end;
    clamp_hits(k) = out.clamp_hits;

    m_nom = out.metr;
    PFA_nom(k)    = m_nom.PFA_top;
    IDR_nom(k)    = m_nom.IDR_max;
    dP95_nom(k)   = m_nom.dP_orf_q95;
    Qcap95_nom(k) = m_nom.Qcap_ratio_q95;
    cav_nom(k)    = m_nom.cav_pct;

    m_w = out.weighted;
    PFA_w(k)    = m_w.PFA_top;
    IDR_w(k)    = m_w.IDR_max;
    dP95_w(k)   = m_w.dP_orf_q95;
    Qcap95_w(k) = m_w.Qcap_ratio_q95;
    if isfield(m_w,'Q_q95'),        Q_q95_w(k)  = m_w.Q_q95; end
    if isfield(m_w,'Q_q50'),        Q_q50_w(k)  = m_w.Q_q50; end
    if isfield(m_w,'dP_orf_q50'),   dP50_w(k)   = m_w.dP_orf_q50; end

    m_ws = out.worst;
    PFA_worst(k)    = m_ws.PFA_top;
    IDR_worst(k)    = m_ws.IDR_max;
    dP95_worst(k)   = m_ws.dP_orf_q95;
    Qcap95_worst(k) = m_ws.Qcap_ratio_q95;
    if isfield(m_ws,'Q_q95'),           Q_q95_worst(k)        = m_ws.Q_q95; end
    if isfield(m_ws,'Q_q50'),           Q_q50_worst(k)        = m_ws.Q_q50; end
    if isfield(m_ws,'dP_orf_q50'),      dP_orf_q50_worst(k)   = m_ws.dP_orf_q50; end
    if isfield(m_ws,'PF_p95'),          PF_p95_worst(k)       = m_ws.PF_p95; end
    which_mu_PFA(k) = m_ws.which_mu.PFA_top;
    which_mu_IDR(k) = m_ws.which_mu.IDR_max;
    T_end_worst(k)  = m_ws.T_oil_end;
    mu_end_worst(k) = m_ws.mu_end;
    cav_pct_worst(k)= m_ws.cav_pct;
    % New: assign worst-case top story envelope metrics
    if isfield(m_ws,'x10_max_D'),        x10_max_D_worst(k)    = m_ws.x10_max_D; end
    if isfield(m_ws,'a10abs_max_D'),     a10abs_max_D_worst(k) = m_ws.a10abs_max_D; end
    % Fallbacks from single-window metrics if worst lacks the renamed fields
    if isfield(m_ws,'x10_pk_D') && ~isfield(m_ws,'x10_max_D')
        x10_max_D_worst(k) = m_ws.x10_pk_D;
    end
    if isfield(m_ws,'a10abs_pk_D') && ~isfield(m_ws,'a10abs_max_D')
        a10abs_max_D_worst(k) = m_ws.a10abs_pk_D;
    end
    if isfield(m_ws,'E_orifice_full'),    E_orifice_sum(k) = m_ws.E_orifice_full; end
    if isfield(m_ws,'E_struct_full'),     E_struct_sum(k)  = m_ws.E_struct_full; end
    if isfield(m_ws,'E_ratio_full'),      E_ratio(k)       = m_ws.E_ratio_full; end
    qc_all_mu(k)    = out.qc_all_mu;

    valsPFA = arrayfun(@(s) s.metr.PFA_top, out.mu_results);
    [maxPFA_rec, idxPFA] = max(valsPFA);
    if maxPFA_rec > worstPFA
        worstPFA = maxPFA_rec;
        worstPFA_name = out.name;
        worstPFA_mu   = out.mu_results(idxPFA).mu_factor;
    end
    valsIDR = arrayfun(@(s) s.metr.IDR_max, out.mu_results);
    [maxIDR_rec, idxIDR] = max(valsIDR);
    if maxIDR_rec > worstIDR
        worstIDR = maxIDR_rec;
        worstIDR_name = out.name;
        worstIDR_mu   = out.mu_results(idxIDR).mu_factor;
    end
end

summary = struct();

summary.table = table(names, scale, SaT1, t5, t95, coverage, rank_score, policy_col, order_col, cooldown_col, ...
    PFA_nom, IDR_nom, dP95_nom, Qcap95_nom, cav_nom, ...
    PFA_w, IDR_w, dP95_w, Qcap95_w, Q_q95_w, Q_q50_w, dP50_w, ...
    PFA_worst, IDR_worst, dP95_worst, dP_orf_q50_worst, Qcap95_worst, Q_q95_worst, Q_q50_worst, PF_p95_worst, ...
    cav_pct_worst, x10_max_D_worst, a10abs_max_D_worst, E_orifice_sum, E_struct_sum, E_ratio, ...
    which_mu_PFA, which_mu_IDR, T_end_worst, mu_end_worst, qc_all_mu, ...
    T_start, T_end, mu_end, clamp_hits, ...
    'VariableNames', {'name','scale','SaT1','t5','t95','coverage','rank_score','policy','order','cooldown_s', ...
    'PFA_nom','IDR_nom','dP95_nom','Qcap95_nom','cav_nom', ...
    'PFA_w','IDR_w','dP95_w','Qcap95_w','Q_q95_w','Q_q50_w','dP50_w', ...
    'PFA_worst','IDR_worst','dP95_worst','dP_orf_q50_worst','Qcap95_worst','Q_q95_worst','Q_q50_worst','PF_p95_worst', ...
    'cav_pct_worst','x10_max_D_worst','a10abs_max_D_worst','E_orifice_sum','E_struct_sum','E_ratio', ...
    'which_mu_PFA','which_mu_IDR','T_end_worst','mu_end_worst','qc_all_mu', ...
    'T_start','T_end','mu_end','clamp_hits'});

% Aliases for compatibility with consumers
summary.table.T_oil_end_worst = summary.table.T_end_worst;
summary.table.dP_orf_q95_worst = summary.table.dP95_worst;
try
    summary.table.energy_tot_sum = summary.table.E_orifice_sum + summary.table.E_struct_sum;
catch
end
summary.all_out = all_out;

% --- QC flags and reason codes for summary.csv consumers ---
% Use thresholds from opts if provided, else defaults consistent with runners
thr_default = struct('dP95_max',50e6,'Qcap95_max',0.5,'cav_pct_max',0,'T_end_max',75,'mu_end_min',0.5);
if isfield(opts,'thr') && ~isempty(opts.thr)
    thr = opts.thr;
    fns = fieldnames(thr_default);
    for ii=1:numel(fns)
        if ~isfield(thr,fns{ii}) || isempty(thr.(fns{ii}))
            thr.(fns{ii}) = thr_default.(fns{ii});
        end
    end
else
    thr = thr_default;
end
ok_T    = summary.table.T_end_worst   <= thr.T_end_max;
ok_mu   = summary.table.mu_end_worst  >= thr.mu_end_min;
ok_dP   = summary.table.dP95_worst    <= thr.dP95_max;
ok_Qcap = summary.table.Qcap95_worst  <  thr.Qcap95_max;
ok_cav  = summary.table.cav_pct_worst == 0;
qc_reason = strings(height(summary.table),1);
for r = 1:height(summary.table)
    bad = {};
    if ~ok_T(r),    bad{end+1}='T';  end %#ok<AGROW>
    if ~ok_mu(r),   bad{end+1}='mu'; end %#ok<AGROW>
    if ~ok_dP(r),   bad{end+1}='dP'; end %#ok<AGROW>
    if ~ok_Qcap(r), bad{end+1}='Qcap'; end %#ok<AGROW>
    if ~ok_cav(r),  bad{end+1}='cav'; end %#ok<AGROW>
    qc_reason(r) = strjoin(bad,',');
end
summary.table.ok_T = ok_T;
summary.table.ok_mu = ok_mu;
summary.table.ok_dP = ok_dP;
summary.table.ok_Qcap = ok_Qcap;
summary.table.ok_cav = ok_cav;
summary.table.qc_reason = qc_reason;

if ~isfield(opts,'quiet') || ~opts.quiet
    fprintf('Worst PFA: %s, mu=%.2f\n', worstPFA_name, worstPFA_mu);
    fprintf('Worst IDR: %s, mu=%.2f\n', worstIDR_name, worstIDR_mu);
end

if do_export
    export_results(outdir, scaled, params, opts, summary, all_out);
    if ~Utils.getfield_default(opts,'quiet',false)
        diary off;
    end
end

end


function P = run_policies_windowed(scaled, params, opts)
%RUN_POLICIES_WINDOWED Evaluate thermal reset policies and record orders.
%   P = RUN_POLICIES_WINDOWED(SCALED, PARAMS, OPTS) orchestrates calls to
%   RUN_BATCH_WINDOWED for combinations of thermal reset policies and record
%   processing orders.  OPTS.policies selects a subset of {'each','carry',
%   'cooldown'} and OPTS.orders selects from {'natural','random','worst_first'}.
%   When 'cooldown' is included, OPTS.cooldown_s_list specifies the cooldown
%   durations to test.  The resulting struct array P contains, for each
%   combination, the policy, order, cooldown duration, summary table, QC stats
%   and deviations relative to the baseline each/natural run.
%
%   OPTS.mu_factors and OPTS.mu_weights mirror RUN_BATCH_WINDOWED defaults.
%   OPTS.rng_seed controls reproducibility of the 'random' order.
%
%   Example:
%       opts.policies = {'each','carry','cooldown'};
%       opts.orders   = {'natural','worst_first'};
%       opts.cooldown_s_list = [60 180 300];
%       P = run_policies_windowed(scaled, params, opts);
%
%   See RUN_BATCH_WINDOWED for additional options.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'mu_factors'), opts.mu_factors = [0.75 1.00 1.25]; end
if ~isfield(opts,'mu_weights'), opts.mu_weights = [0.2 0.6 0.2]; end
if ~isfield(opts,'policies'), opts.policies = {'each','carry','cooldown'}; end
if ~isfield(opts,'orders'), opts.orders = {'natural','random','worst_first'}; end
if ~isfield(opts,'cooldown_s_list'), opts.cooldown_s_list = [60 180 300]; end
if ~isfield(opts,'rng_seed'), opts.rng_seed = 42; end
if ~isfield(opts,'rank_metric'), opts.rank_metric = 'E_orifice_win'; end

% QC thresholds for logging (kept in opts for downstream use)
thr_default = struct('dP95_max',50e6,'Qcap95_max',0.5,'cav_pct_max',0,'T_end_max',75,'mu_end_min',0.5);
if ~isfield(opts,'thr') || isempty(opts.thr)
    opts.thr = thr_default;
else
    fns = fieldnames(thr_default);
    for ii=1:numel(fns)
        if ~isfield(opts.thr,fns{ii}) || isempty(opts.thr.(fns{ii}))
            opts.thr.(fns{ii}) = thr_default.(fns{ii});
        end
    end
end

% Quiet/export flags and output dir
quiet = isfield(opts,'quiet') && opts.quiet;
% Default to exporting results for policy runs unless explicitly disabled
do_export = ~isfield(opts,'do_export') || opts.do_export;
ts = datestr(now,'yyyymmdd_HHMMSS_FFF'); outdir = fullfile('out', ts);
if do_export && ~quiet
    if ~exist(outdir,'dir'), mkdir(outdir); end
    diary(fullfile(outdir,'console.log'));
else
    if ~exist(outdir,'dir'), mkdir(outdir); end
end

% Console header
try
    imode = Utils.getfield_default(opts,'IM_mode','');
    band  = Utils.getfield_default(opts,'band_fac',[NaN NaN]);
    sb    = Utils.getfield_default(opts,'s_bounds',[NaN NaN]);
    if ~quiet
        fprintf('Run @ %s | outdir=%s\n', ts, outdir);
        fprintf('IM_mode=%s, band=[%.3g,%.3g], s_bounds=[%.2f,%.2f]\n', imode, band(1), band(2), sb(1), sb(2));
        fprintf('policies=%s | orders=%s | cooldown_s_list=%s | rng_seed=%d\n', ...
            strjoin(opts.policies,','), strjoin(opts.orders,','), sprintf('%d ', opts.cooldown_s_list), opts.rng_seed);
        if isfield(opts,'TRIM_names') && ~isempty(opts.TRIM_names)
            fprintf('TRIM: %s\n', strjoin(opts.TRIM_names,', '));
        end
        fprintf('QC thr: dP95<=%.1f MPa, Qcap95<%.2f, cav%%=%g, T_end<=%g C, mu_end>=%0.2f Pa*s\n', ...
            opts.thr.dP95_max/1e6, opts.thr.Qcap95_max, opts.thr.cav_pct_max*100, opts.thr.T_end_max, opts.thr.mu_end_min);
    end
catch ME
    warning('run_batch_windowed header: %s', ME.message);
end

% Echo hydraulic/thermal key params
try
    n_orf = NaN; if isfield(params,'A_o'), n_orf = numel(params.A_o); end
    A_o = Utils.getfield_default(params,'A_o',NaN);
    d_o = NaN; try, d_o = sqrt(4*mean(A_o)/pi); catch, end
    Qcap_big = Utils.getfield_default(params,'Qcap_big',NaN);
    hA = NaN; if isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), hA = params.thermal.hA_W_perK; end
    resFactor = Utils.getfield_default(params,'resFactor',NaN);
    tg = Utils.getfield_default(params,'toggle_gain',NaN);
    tgv = tg(:); tgmin = min(tgv); tgmed = median(tgv); tgmax = max(tgv);
    if ~quiet
        fprintf('Hydraulics: n_orf=%g, d_o~=%g m, A_o=%s, Qcap_big=%g, hA=%g, resFactor=%g, toggle_gain[min/med/max]=[%g %g %g]\n', ...
            n_orf, d_o, mat2str(size(A_o)), Qcap_big, hA, resFactor, tgmin, tgmed, tgmax);
    end
catch ME
    warning('run_batch_windowed echo params: %s', ME.message);
end

nRec = numel(scaled);

% Baseline run for deltas and worst_first ordering
base_opts = opts; base_opts.thermal_reset = 'each'; base_opts.order = 'natural';
base_opts.do_export = false;
[base_summary, base_all] = run_batch_windowed(scaled, params, base_opts);
% Baseline aggregates and log
basePFA_w_mean = mean(base_summary.table.PFA_w);
baseIDR_w_mean = mean(base_summary.table.IDR_w);
baseTend_worst = max(base_summary.table.T_end_worst);
base_dP95_worst_max = max(base_summary.table.dP95_worst);
base_Qcap95_worst_max = max(base_summary.table.Qcap95_worst);
base_cav_worst_max = max(base_summary.table.cav_pct_worst);
base_mu_end_worst_min = min(base_summary.table.mu_end_worst);
base_qc_pass = sum(base_summary.table.qc_all_mu);
base_qc_n = height(base_summary.table);
if ~quiet
    fprintf(['BASE (each/natural): PFA_w=%.3g, IDR_w=%.3g, dP95_worst=%.3g MPa, Qcap95_worst=%.2f, ' ...
            'cav%%_worst=%.1f, T_end_worst=%.1f C, mu_end_worst=%.2f, qc_rate=%d/%d\n'], ...
        basePFA_w_mean, baseIDR_w_mean, base_dP95_worst_max/1e6, base_Qcap95_worst_max, base_cav_worst_max*100, baseTend_worst, base_mu_end_worst_min, base_qc_pass, base_qc_n);
end

% Pre-compute orders
orders_struct.natural = 1:nRec;
if any(strcmp(opts.orders,'random'))
    rng(opts.rng_seed);
    orders_struct.random = randperm(nRec);
end
if any(strcmp(opts.orders,'worst_first'))
    switch lower(opts.rank_metric)
        case 'pfa_top'
            rk = cellfun(@(s) s.metr.PFA_top, base_all);
        case 'idr_max'
            rk = cellfun(@(s) s.metr.IDR_max, base_all);
        case 'pfa_w'
            rk = cellfun(@(s) s.weighted.PFA_top, base_all);
        otherwise
            rk = cellfun(@(s) s.metr.E_orifice_win, base_all);
    end
    [~,idx] = sort(rk,'descend');
    orders_struct.worst_first = idx;
    try
        rank_names = {scaled(idx).name};
        if ~quiet
            fprintf('worst_first ranking by %s: %s\n', opts.rank_metric, strjoin(rank_names, ', '));
        end
    catch ME
        warning('run_batch_windowed ranking: %s', ME.message);
    end
end

% Iterate combinations
P = struct('policy',{},'order',{},'cooldown_s',{},'summary',{},'qc',{},'deltas',{});
for ip = 1:numel(opts.policies)
    pol = opts.policies{ip};
    for io = 1:numel(opts.orders)
        ord = opts.orders{io};
        if strcmp(pol,'cooldown')
            cds = opts.cooldown_s_list(:)';
        else
            cds = NaN;
        end
        for ic = 1:numel(cds)
            cdval = cds(ic);
            perm = orders_struct.(ord);
            scaled_run = scaled(perm);
            run_opts = opts;
            if isfield(run_opts,'cooldown_s'), run_opts = rmfield(run_opts,'cooldown_s'); end
            run_opts.order = ord;
            run_opts.thermal_reset = pol;
            run_opts.do_export = false;
            if strcmp(pol,'cooldown'), run_opts.cooldown_s = cdval; end
            [summary, ~] = run_batch_windowed(scaled_run, params, run_opts);

            qc.pass_fraction = mean(summary.table.qc_all_mu);
            qc.n = height(summary.table);

            % Aggregates for this combination
            curPFA_w_mean = mean(summary.table.PFA_w);
            curIDR_w_mean = mean(summary.table.IDR_w);
            curTend_worst = max(summary.table.T_end_worst);
            deltas = struct('PFA_w', curPFA_w_mean - basePFA_w_mean, ...
                            'IDR_w', curIDR_w_mean - baseIDR_w_mean, ...
                            'T_end_worst', curTend_worst - baseTend_worst);

            % log worst cases
            [worstPFA, idxP] = max(summary.table.PFA_worst);
            nP = summary.table.name{idxP};
            muP = summary.table.which_mu_PFA(idxP);
            if ~quiet, fprintf('Worst PFA (%s,%s,mu=%.2f): %s\n', pol, ord, muP, nP); end
            [worstIDR, idxI] = max(summary.table.IDR_worst); %#ok<NASGU>
            nI = summary.table.name{idxI};
            muI = summary.table.which_mu_IDR(idxI);
            if ~quiet, fprintf('Worst IDR (%s,%s,mu=%.2f): %s\n', pol, ord, muI, nI); end

            % policy comparison rule logging against baseline each/natural
            tolPFA = 0.15 * basePFA_w_mean;
            tolIDR = 0.15 * baseIDR_w_mean;
            passPFA = abs(deltas.PFA_w) <= tolPFA;
            passIDR = abs(deltas.IDR_w) <= tolIDR;
            if ~quiet
                fprintf('Delta vs base: dPFA_w=%.4g (|d|<=%.4g? %d), dIDR_w=%.4g (|d|<=%.4g? %d), qc_rate=%.2f\n', ...
                    deltas.PFA_w, tolPFA, passPFA, deltas.IDR_w, tolIDR, passIDR, qc.pass_fraction);
            end

            % Single-line summary for this combination
            n_pass = sum(summary.table.qc_all_mu);
            n_tot  = height(summary.table);
            PFAw   = curPFA_w_mean; IDRw = curIDR_w_mean; T_end_worst_max = curTend_worst;
            pass_flag_15pct = 'OK'; if ~(passPFA && passIDR), pass_flag_15pct = 'FAIL'; end
            if ~quiet
                fprintf(['policy=%s | order=%s | cd=%ds | PFA_w=%.3g (d=%+.2f%%) | IDR_w=%.3g (d=%+.2f%%) | ' ...
                        'T_end_worst=%.1f C | qc_rate=%d/%d %s\n'], ...
                    pol, ord, cdval, PFAw, 100*(PFAw-basePFA_w_mean)/max(basePFA_w_mean,eps), ...
                    IDRw, 100*(IDRw-baseIDR_w_mean)/max(baseIDR_w_mean,eps), T_end_worst_max, n_pass, n_tot, pass_flag_15pct);
                % Echo worst two records by PFA_worst
                try
                    kshow = min(2, height(summary.table));
                    [~,ix] = maxk(summary.table.PFA_worst, kshow);
                    if kshow==2
                        fprintf('  worst2(PFA): %s | %s\n', summary.table.name{ix(1)}, summary.table.name{ix(2)});
                    elseif kshow==1
                        fprintf('  worst2(PFA): %s | -\n', summary.table.name{ix(1)});
                    end
                catch ME
                    warning('run_batch_windowed worst2: %s', ME.message);
                end
            end

            % Optional clamp summary
            try
                if ~quiet
                    tot_clamps = sum(summary.table.clamp_hits);
                    nz_idx = find(summary.table.clamp_hits>0);
                    nz_names = summary.table.name(nz_idx);
                    if ~isempty(nz_names)
                        fprintf('clamp_hits total=%d | records: %s\n', tot_clamps, strjoin(nz_names.', ', '));
                    else
                        fprintf('clamp_hits total=%d\n', tot_clamps);
                    end
                end
            catch ME
                warning('run_batch_windowed clamp summary: %s', ME.message);
            end

            P(end+1) = struct('policy',pol,'order',ord,'cooldown_s',cdval, ...
                'summary',summary.table,'qc',qc,'deltas',deltas); %#ok<AGROW>
        end
    end
end

% Export and close
if do_export
    export_results(outdir, scaled, params, opts, base_summary, base_all, P);
end
if ~quiet, fprintf('Saved to %s\n', outdir); end
if do_export && ~quiet, diary off; end
end
