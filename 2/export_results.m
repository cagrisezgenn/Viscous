function export_results(outdir, scaled, params, opts, summary, all_out, varargin)
%EXPORT_RESULTS Save run outputs to disk.
%   EXPORT_RESULTS(OUTDIR,SCALED,PARAMS,OPTS,SUMMARY,ALL_OUT) stores
%   snapshot, summary tables and per-record details under OUTDIR.  Optional
%   policy results can be supplied as an extra argument.

if nargin < 7, varargin = {}; end

%% Meta / snapshot
try
    save(fullfile(outdir,'snapshot.mat'), 'params','opts','scaled','-v7.3');
    idx = struct('name',{scaled.name}, 'dt',{scaled.dt}, 'dur',{scaled.dur}, ...
                 'PGA',{scaled.PGA}, 'PGV',{scaled.PGV}, 'IM',{scaled.IM}, ...
                 'scale',{scaled.scale});
    writetable(struct2table(idx), fullfile(outdir,'scaled_index.csv'));
catch ME
    warning('export_results:snapshot', ME.message);
end

%% Summaries
try
    writetable(summary.table, fullfile(outdir,'summary.csv'));
catch ME
    warning('export_results:summary', ME.message);
end
try
    save(fullfile(outdir,'summary_full.mat'), 'summary','all_out','-v7.3');
catch ME
    warning('export_results:summary_full', ME.message);
end
try
    pass = sum(summary.table.qc_all_mu);
    fail = height(summary.table) - pass;
    qc_counts = table(pass, fail);
    [~, idx] = maxk(summary.table.PFA_worst, min(3,height(summary.table)));
    worst_tbl = summary.table(idx, {'name','PFA_worst','IDR_worst'});
    qc_file = fullfile(outdir,'qc_summary.csv');
    writetable(qc_counts, qc_file);
    writetable(worst_tbl, qc_file, 'WriteMode','append');
catch ME
    warning('export_results:qc', ME.message);
end

%% Record-based details
for k = 1:numel(all_out)
    out = all_out{k};
    recdir = fullfile(outdir, ['rec-' sanitize_name(out.name)]);
    if ~exist(recdir,'dir'), mkdir(recdir); end
    % window
    try
        w = struct('t5', out.win.t5, 't95', out.win.t95, ...
                   'pad', getfield(out.win,'pad',[]), ...
                   'coverage', out.win.coverage);
        writejson(fullfile(recdir,'window.json'), w);
    catch ME
        warning('export_results:window', ME.message);
    end
    % mu results
    try
        if isfield(out, 'mu_results') && ~isempty(out.mu_results)
            mu_tbl = table([out.mu_results.mu_factor]', ...
                           arrayfun(@(s) s.metr.PFA_top, out.mu_results)', ...
                           arrayfun(@(s) s.metr.IDR_max, out.mu_results)', ...
                           arrayfun(@(s) s.metr.dP_orf_q95, out.mu_results)', ...
                           arrayfun(@(s) s.metr.Qcap_ratio_q95, out.mu_results)', ...
                           arrayfun(@(s) s.metr.T_oil_end, out.mu_results)', ...
                           arrayfun(@(s) s.metr.mu_end, out.mu_results)', ...
                           'VariableNames', {'mu','PFA_top','IDR_max','dP95','Qcap95','T_end','mu_end'});
            writetable(mu_tbl, fullfile(recdir,'mu_results.csv'));
        end
    catch ME
        warning('export_results:mu', ME.message);
    end
    % nominal metrics
    try
        m = out.metr;
        m_tbl = table(m.PFA_top, m.IDR_max, m.dP_orf_q95, m.Qcap_ratio_q95, ...
                      m.cav_pct, m.E_orifice_win, ...
                      getfield(m,'E_struct_win',NaN), ...
                      getfield(m,'E_mech_win',NaN), ...
                      out.T_end, out.mu_end, ...
                      'VariableNames', {'PFA_top','IDR_max','dP95','Qcap95','cav_pct', ...
                      'E_orifice','E_struct','E_mech','T_end','mu_end'});
        writetable(m_tbl, fullfile(recdir,'metrics_win.csv'));
    catch ME
        warning('export_results:metrics', ME.message);
    end
    % time series downsample
    try
        if isfield(out,'ts') && ~isempty(out.ts)
            ds = 5;
            if isfield(opts,'export') && isfield(opts.export,'ds')
                ds = opts.export.ds;
            end
            ts_ds = downsample_ts(out.ts, ds);
            save(fullfile(recdir,'ts_ds.mat'),'ts_ds','-v7.3');
        end
    catch ME
        warning('export_results:ts', ME.message);
    end
    % optional plots
    if isfield(opts,'export') && isfield(opts.export,'plots') && opts.export.plots
        try
            figdir = fullfile(recdir,'figures');
            if ~exist(figdir,'dir'), mkdir(figdir); end
            if isfield(out,'mu_results') && ~isempty(out.mu_results)
                muvals = [out.mu_results.mu_factor]';
                PFAvals = arrayfun(@(s) s.metr.PFA_top, out.mu_results)';
                IDRvals = arrayfun(@(s) s.metr.IDR_max, out.mu_results)';
                f = figure('visible','off');
                bar(muvals, PFAvals);
                xlabel('\mu'); ylabel('PFA');
                saveas(f, fullfile(figdir,'PFA.png'));
                close(f);
                f = figure('visible','off');
                bar(muvals, IDRvals);
                xlabel('\mu'); ylabel('IDR');
                saveas(f, fullfile(figdir,'IDR.png'));
                close(f);
            end
        catch ME
            warning('export_results:plots', ME.message);
        end
    end
end

%% Policy results
if ~isempty(varargin)
    P = varargin{1};
    try
        idx_tbl = table();
        for i = 1:numel(P)
            pol = P(i);
            fname = sprintf('policy_%s_%s_cd%s.csv', ...
                sanitize_name(pol.policy), sanitize_name(pol.order), num2str(pol.cooldown_s));
            writetable(pol.summary, fullfile(outdir, fname));
            idx_tbl = [idx_tbl; table({pol.policy},{pol.order},pol.cooldown_s, ...
                pol.qc.pass_fraction, max(pol.summary.PFA_w), max(pol.summary.IDR_w), ...
                max(pol.summary.PFA_worst), max(pol.summary.IDR_worst), ...
                'VariableNames', {'policy','order','cooldown_s','qc_pass_fraction','PFA_w','IDR_w','PFA_worst','IDR_worst'})]; %#ok<AGROW>
        end
        writetable(idx_tbl, fullfile(outdir,'policy_index.csv'));
    catch ME
        warning('export_results:policy', ME.message);
    end
end
end

