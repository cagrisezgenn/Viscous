function export_results(outdir, scaled, params, opts, summary, all_out, varargin)
%EXPORT_RESULTS Persist results to out/<timestamp>/ directory.
%   EXPORT_RESULTS(OUTDIR, SCALED, PARAMS, OPTS, SUMMARY, ALL_OUT)
%   writes various logs and per-record details.  Optional additional
%   structures (e.g., policy results) may be passed in VARARGIN.

if nargin < 1 || isempty(outdir), return; end
if ~exist(outdir,'dir'), mkdir(outdir); end

%% Meta/snapshot
try
    save(fullfile(outdir,'snapshot.mat'), 'params','opts','scaled','-v7.3');
catch
end
try
    names = {scaled.name}';
    dt = arrayfun(@(s) getfield_default(s,'dt',NaN), scaled)';
    dur = arrayfun(@(s) getfield_default(s,'duration',NaN), scaled)';
    PGA = arrayfun(@(s) getfield_default(s,'PGA',NaN), scaled)';
    PGV = arrayfun(@(s) getfield_default(s,'PGV',NaN), scaled)';
    IM  = arrayfun(@(s) getfield_default(s,'IM',NaN), scaled)';
    sc  = arrayfun(@(s) getfield_default(s,'scale',NaN), scaled)';
    tbl = table(names, dt, dur, PGA, PGV, IM, sc, ...
        'VariableNames',{'name','dt','dur','PGA','PGV','IM','scale'});
    writetable(tbl, fullfile(outdir,'scaled_index.csv'));
catch
end

%% Summaries
try
    writetable(summary.table, fullfile(outdir,'summary.csv'));
catch
end
try
    save(fullfile(outdir,'summary_full.mat'), 'summary','all_out','-v7.3');
catch
end
try
    pass = sum(summary.table.qc_all_mu);
    fail = height(summary.table) - pass;
    [~,idx] = maxk(summary.table.PFA_worst, min(3,height(summary.table)));
    worst = summary.table.name(idx);
    qc_tbl = table(pass, fail, worst, 'VariableNames',{'pass','fail','worst'});
    writetable(qc_tbl, fullfile(outdir,'qc_summary.csv'));
catch
end

%% Record-based details
for k = 1:numel(all_out)
    out = all_out{k};
    try
        recdir = fullfile(outdir, ['rec-' sanitize_name(out.name)]);
        if ~exist(recdir,'dir'), mkdir(recdir); end
        % window.json
        try
            w = out.win;
            win_struct = struct('t5',w.t5,'t95',w.t95,'pad',getfield_default(w,'pad',0), ...
                                'coverage',w.coverage);
            writejson(win_struct, fullfile(recdir,'window.json'));
        catch
        end
        % mu_results.csv
        if isfield(out,'mu_results') && ~isempty(out.mu_results)
            try
                mu_tbl = struct2table(arrayfun(@(s) s.metr, out.mu_results));
                mu_tbl.mu = [out.mu_results.mu_factor]';
                mu_tbl = movevars(mu_tbl,'mu','before',1);
                writetable(mu_tbl, fullfile(recdir,'mu_results.csv'));
            catch
            end
        end
        % metrics_win.csv
        try
            m = out.metr;
            mstruct = struct('PFA_top',m.PFA_top,'IDR_max',m.IDR_max, ...
                             'dP95',m.dP_orf_q95,'Qcap95',m.Qcap_ratio_q95, ...
                             'cav_pct',m.cav_pct,'T_end',out.T_end,'mu_end',out.mu_end);
            if isfield(m,'E_orifice_win'), mstruct.E_orf_win = m.E_orifice_win; end
            if isfield(m,'E_struct_win'), mstruct.E_struct_win = m.E_struct_win; end
            metr_tbl = struct2table(mstruct);
            writetable(metr_tbl, fullfile(recdir,'metrics_win.csv'));
        catch
        end
        % ts_ds.mat
        if isfield(out,'ts') && ~isempty(out.ts)
            try
                ds = 5;
                if isfield(opts,'export') && isfield(opts.export,'ds')
                    ds = opts.export.ds;
                end
                ts_ds = downsample_ts(out.ts, ds);
                save(fullfile(recdir,'ts_ds.mat'),'ts_ds','-v7.3');
            catch
            end
        end
        % optional plots
        if isfield(opts,'export') && isfield(opts.export,'plots') && opts.export.plots
            try
                figdir = fullfile(recdir,'figures');
                if ~exist(figdir,'dir'), mkdir(figdir); end
                f = figure('Visible','off');
                bar(out.metr.IDR_max); title('IDR_{max}');
                saveas(f, fullfile(figdir,'IDR_max.png')); close(f);
                f = figure('Visible','off');
                bar(out.metr.PFA_top); title('PFA_{top}');
                saveas(f, fullfile(figdir,'PFA_top.png')); close(f);
            catch
            end
        end
    catch
    end
end

%% Policy results
if ~isempty(varargin)
    P = varargin{1};
    try
        pol = {P.policy}';
        ord = {P.order}';
        cd  = [P.cooldown_s]';
        qc_pass = arrayfun(@(s) s.qc.pass_fraction, P)';
        qc_n    = arrayfun(@(s) s.qc.n, P)';
        PFA_w   = arrayfun(@(s) max(s.summary.PFA_nom), P)';
        IDR_w   = arrayfun(@(s) max(s.summary.IDR_nom), P)';
        PFA_ws  = arrayfun(@(s) max(s.summary.PFA_worst), P)';
        IDR_ws  = arrayfun(@(s) max(s.summary.IDR_worst), P)';
        idx_tbl = table(pol, ord, cd, qc_pass, qc_n, PFA_w, IDR_w, PFA_ws, IDR_ws, ...
            'VariableNames',{'policy','order','cooldown_s','qc_pass_frac','qc_n', ...
                             'PFA_w','IDR_w','PFA_worst','IDR_worst'});
        writetable(idx_tbl, fullfile(outdir,'policy_index.csv'));
        for i = 1:numel(P)
            fname = sprintf('policy_%s_%s_cd%s.csv', ...
                sanitize_name(P(i).policy), sanitize_name(P(i).order), ...
                sanitize_name(num2str(P(i).cooldown_s)));
            writetable(P(i).summary, fullfile(outdir,fname));
        end
    catch
    end
end
end
