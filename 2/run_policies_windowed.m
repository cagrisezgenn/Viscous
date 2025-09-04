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

% Always create outdir and start diary for this run
ts = datestr(now,'yyyymmdd_HHMMSS'); outdir = fullfile('out', ts);
if ~exist(outdir,'dir'), mkdir(outdir); end
diary(fullfile(outdir,'console.log'));

% Default to exporting results for policy runs unless explicitly disabled
do_export = ~isfield(opts,'do_export') || opts.do_export;

% Console header
try
    imode = getfield_default(opts,'IM_mode','');
    band  = getfield_default(opts,'band_fac',[NaN NaN]);
    sb    = getfield_default(opts,'s_bounds',[NaN NaN]);
    fprintf('Run @ %s | outdir=%s\n', ts, outdir);
    fprintf('IM_mode=%s, band=[%.3g,%.3g], s_bounds=[%.2f,%.2f]\n', imode, band(1), band(2), sb(1), sb(2));
    fprintf('policies=%s | orders=%s | cooldown_s_list=%s | rng_seed=%d\n', ...
        strjoin(opts.policies,','), strjoin(opts.orders,','), sprintf('%d ', opts.cooldown_s_list), opts.rng_seed);
    if isfield(opts,'TRIM_names') && ~isempty(opts.TRIM_names)
        fprintf('TRIM: %s\n', strjoin(opts.TRIM_names,', '));
    end
    fprintf('QC thr: dP95<=%.1f MPa, Qcap95<%.2f, cav%%=%g, T_end<=%g C, mu_end>=%0.2f Pa*s\n', ...
        opts.thr.dP95_max/1e6, opts.thr.Qcap95_max, opts.thr.cav_pct_max*100, opts.thr.T_end_max, opts.thr.mu_end_min);
catch
end

% Echo hydraulic/thermal key params
try
    n_orf = NaN; if isfield(params,'A_o'), n_orf = numel(params.A_o); end
    A_o = getfield_default(params,'A_o',NaN);
    d_o = NaN; try, d_o = sqrt(4*mean(A_o)/pi); catch, end
    Qcap_big = getfield_default(params,'Qcap_big',NaN);
    hA = NaN; if isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), hA = params.thermal.hA_W_perK; end
    resFactor = getfield_default(params,'resFactor',NaN);
    tg = getfield_default(params,'toggle_gain',NaN);
    tgv = tg(:); tgmin = min(tgv); tgmed = median(tgv); tgmax = max(tgv);
    fprintf('Hydraulics: n_orf=%g, d_o~=%g m, A_o=%s, Qcap_big=%g, hA=%g, resFactor=%g, toggle_gain[min/med/max]=[%g %g %g]\n', ...
        n_orf, d_o, mat2str(size(A_o)), Qcap_big, hA, resFactor, tgmin, tgmed, tgmax);
catch
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
fprintf(['BASE (each/natural): PFA_w=%.3g, IDR_w=%.3g, dP95_worst=%.3g MPa, Qcap95_worst=%.2f, ' ...
        'cav%%_worst=%.1f, T_end_worst=%.1f C, mu_end_worst=%.2f, qc_rate=%d/%d\n'], ...
    basePFA_w_mean, baseIDR_w_mean, base_dP95_worst_max/1e6, base_Qcap95_worst_max, base_cav_worst_max*100, baseTend_worst, base_mu_end_worst_min, base_qc_pass, base_qc_n);

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
        fprintf('worst_first ranking by %s: %s\n', opts.rank_metric, strjoin(rank_names, ', '));
    catch, end
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
            fprintf('Worst PFA (%s,%s,mu=%.2f): %s\n', pol, ord, muP, nP);
            [worstIDR, idxI] = max(summary.table.IDR_worst); %#ok<NASGU>
            nI = summary.table.name{idxI};
            muI = summary.table.which_mu_IDR(idxI);
            fprintf('Worst IDR (%s,%s,mu=%.2f): %s\n', pol, ord, muI, nI);

            % policy comparison rule logging against baseline each/natural
            tolPFA = 0.15 * basePFA_w_mean;
            tolIDR = 0.15 * baseIDR_w_mean;
            passPFA = abs(deltas.PFA_w) <= tolPFA;
            passIDR = abs(deltas.IDR_w) <= tolIDR;
            fprintf('Delta vs base: dPFA_w=%.4g (|d|<=%.4g? %d), dIDR_w=%.4g (|d|<=%.4g? %d), qc_rate=%.2f\n', ...
                deltas.PFA_w, tolPFA, passPFA, deltas.IDR_w, tolIDR, passIDR, qc.pass_fraction);

            % Single-line summary for this combination
            n_pass = sum(summary.table.qc_all_mu);
            n_tot  = height(summary.table);
            PFAw   = curPFA_w_mean; IDRw = curIDR_w_mean; T_end_worst_max = curTend_worst;
            pass_flag_15pct = 'OK'; if ~(passPFA && passIDR), pass_flag_15pct = 'FAIL'; end
            fprintf(['policy=%s | order=%s | cd=%ds | PFA_w=%.3g (d=%+.2f%%) | IDR_w=%.3g (d=%+.2f%%) | ' ...
                    'T_end_worst=%.1f C | qc_rate=%d/%d %s\n'], ...
                pol, ord, cdval, PFAw, 100*(PFAw-basePFA_w_mean)/max(basePFA_w_mean,eps), ...
                IDRw, 100*(IDRw-baseIDR_w_mean)/max(baseIDR_w_mean,eps), T_end_worst_max, n_pass, n_tot, pass_flag_15pct);

            % Optional clamp summary
            try
                tot_clamps = sum(summary.table.clamp_hits);
                nz_idx = find(summary.table.clamp_hits>0);
                nz_names = summary.table.name(nz_idx);
                if ~isempty(nz_names)
                    fprintf('clamp_hits total=%d | records: %s\n', tot_clamps, strjoin(nz_names.', ', '));
                else
                    fprintf('clamp_hits total=%d\n', tot_clamps);
                end
            catch
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
fprintf('Saved to %s\n', outdir);
diary off;
end
