function [summary, all_out, P] = run_policies_windowed(scaled, params, opts)
%RUN_POLICIES_WINDOWED Evaluate thermal reset policies and record orders.
%   [SUMMARY, ALL_OUT, P] = RUN_POLICIES_WINDOWED(SCALED, PARAMS, OPTS)
%   orchestrates calls to
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
if ~isfield(opts,'cooldown_s_list'), opts.cooldown_s_list = 60; end
if ~isfield(opts,'rng_seed'), opts.rng_seed = 42; end
if ~isfield(opts,'do_export'), opts.do_export = true; end

if opts.do_export
    ts = datestr(now,'yyyymmdd_HHMMSS');
    outdir = fullfile('out', ts);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    diary(fullfile(outdir,'console.log'));
else
    outdir = '';
end

nRec = numel(scaled);

% Baseline run for deltas and worst_first ordering
base_opts = opts; base_opts.thermal_reset = 'each'; base_opts.order = 'natural';
base_opts.do_export = false;
[base_summary, base_all] = run_batch_windowed(scaled, params, base_opts);
basePFA = max(base_summary.table.PFA_nom);
baseIDR = max(base_summary.table.IDR_nom);
baseTend = max(base_summary.table.T_end);

% Pre-compute orders
orders_struct.natural = 1:nRec;
if any(strcmp(opts.orders,'random'))
    rng(opts.rng_seed);
    orders_struct.random = randperm(nRec);
end
if any(strcmp(opts.orders,'worst_first'))
    E = cellfun(@(s) s.metr.E_orifice_win, base_all);
    [~,idx] = sort(E,'descend');
    orders_struct.worst_first = idx;
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

            curPFA = max(summary.table.PFA_nom);
            curIDR = max(summary.table.IDR_nom);
            curTend = max(summary.table.T_end);
            deltas = struct('PFA', curPFA - basePFA, ...
                            'IDR', curIDR - baseIDR, ...
                            'T_end', curTend - baseTend);

            % log worst cases
            [worstPFA, idxP] = max(summary.table.PFA_worst);
            nP = summary.table.name{idxP};
            muP = summary.table.which_mu_PFA(idxP);
            fprintf('Worst PFA (%s,%s,mu=%.2f): %s\n', pol, ord, muP, nP);
            [worstIDR, idxI] = max(summary.table.IDR_worst); %#ok<NASGU>
            nI = summary.table.name{idxI};
            muI = summary.table.which_mu_IDR(idxI);
            fprintf('Worst IDR (%s,%s,mu=%.2f): %s\n', pol, ord, muI, nI);

            P(end+1) = struct('policy',pol,'order',ord,'cooldown_s',cdval, ...
                'summary',summary.table,'qc',qc,'deltas',deltas); %#ok<AGROW>
        end
    end
end
summary = base_summary;
all_out = base_all;
if opts.do_export
    try
        export_results(outdir, scaled, params, opts, summary, all_out, P);
    catch ME
        warning('export_results failed: %s', ME.message);
    end
    diary off;
end

end
