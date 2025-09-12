function export_results(outdir, scaled, params, opts, summary, all_out, varargin)
%EXPORT_RESULTS Persist results to out/<timestamp>/ directory.
%   EXPORT_RESULTS(OUTDIR, SCALED, PARAMS, OPTS, SUMMARY, ALL_OUT)
%   writes various logs and per-record details.  Optional additional
%   structures (e.g., policy results) may be passed in VARARGIN.

if nargin < 1 || isempty(outdir), return; end
if ~exist(outdir,'dir'), mkdir(outdir); end

%% Meta/snapshot
try
    % IM and trimming info  (BOŞSA BİLE DOLUYA ÇEK!)
    IM_mode  = Utils.getfield_default(opts,'IM_mode','band');
    if isempty(IM_mode), IM_mode = 'band'; end

    band_fac = Utils.getfield_default(opts,'band_fac',[0.8 1.2]);
    if isempty(band_fac), band_fac = [0.8 1.2]; end

    s_bounds = Utils.getfield_default(opts,'s_bounds',[0.5 2.0]);
    if isempty(s_bounds), s_bounds = [0.5 2.0]; end

    TRIM_names = Utils.getfield_default(opts,'TRIM_names',{});
    if isempty(TRIM_names), TRIM_names = {}; end

    % Derived params if available
    params_derived = struct();
    try
        params_derived.C_th = compute_Cth_effective(params);
    catch
    end

    % qc meta (rate + worst3) for snapshot
    try
        pass = sum(summary.table.qc_all_mu);
        fail = height(summary.table) - pass;
        kshow = min(3, height(summary.table));
        [~,ix] = maxk(summary.table.PFA_worst, kshow);
        worst3 = summary.table.name(ix);
        qc_meta = struct('pass',pass,'fail',fail,'worst3',{worst3});
    catch
        qc_meta = struct();
    end

    % snapshot.mat kaydı (P varsa dahil et)
    if ~isempty(varargin)
        P = varargin{1}; %#ok<NASGU>
        save(fullfile(outdir,'snapshot.mat'), ...
             'params','opts','scaled','P', ...
             'IM_mode','band_fac','s_bounds','TRIM_names','params_derived','qc_meta','-v7.3');
    else
        save(fullfile(outdir,'snapshot.mat'), ...
             'params','opts','scaled', ...
             'IM_mode','band_fac','s_bounds','TRIM_names','params_derived','qc_meta','-v7.3');
    end
catch
end

try
    names = {scaled.name}';
    dt = arrayfun(@(s) Utils.getfield_default(s,'dt',NaN), scaled)';
    dur = arrayfun(@(s) Utils.getfield_default(s,'duration',NaN), scaled)';
    PGA = arrayfun(@(s) Utils.getfield_default(s,'PGA',NaN), scaled)';
    PGV = arrayfun(@(s) Utils.getfield_default(s,'PGV',NaN), scaled)';
    IM  = arrayfun(@(s) Utils.getfield_default(s,'IM',NaN), scaled)';
    sc  = arrayfun(@(s) Utils.getfield_default(s,'scale',NaN), scaled)';
    s_cl = arrayfun(@(s) Utils.getfield_default(s,'s_clipped',0), scaled)';
    trimmed = arrayfun(@(s) Utils.getfield_default(s,'trimmed',0), scaled)';
    tbl = table(names, dt, dur, PGA, PGV, IM, sc, s_cl, trimmed, ...
        'VariableNames',{'name','dt','dur','PGA','PGV','IM','scale','s_clipped','trimmed'});
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
        recdir = fullfile(outdir, ['rec-' Utils.sanitize_name(out.name)]);
        if ~exist(recdir,'dir'), mkdir(recdir); end
        % window.json
        try
            w = out.win;
            win_struct = struct('t5',w.t5,'t95',w.t95,'pad',Utils.getfield_default(w,'pad',0), ...
                                'coverage',w.coverage);
            Utils.writejson(win_struct, fullfile(recdir,'window.json'));
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
            mu_mode = Utils.getfield_default(out,'mu_mode','nominal');   % label only
            mu_used = Utils.getfield_default(out,'mu_used',1.0);
            mstruct = struct('PFA_top',m.PFA_top,'IDR_max',m.IDR_max, ...
                             'dP95',m.dP_orf_q95,'Qcap95',m.Qcap_ratio_q95, ...
                             'cav_pct',m.cav_pct,'T_end',out.T_end,'mu_end',out.mu_end);
            % append PF_t_on if present in OUT
            if isfield(out,'PF_t_on')
                mstruct.PF_t_on = out.PF_t_on;
            end
            mstruct.mu_mode = mu_mode; mstruct.mu_used = mu_used;
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
                ts_ds = Utils.downsample_ts(out.ts, ds);
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
        qc_rate = arrayfun(@(s) s.qc.pass_fraction, P)';
        PFA_w_mean = arrayfun(@(s) mean(s.summary.PFA_w), P)';
        IDR_w_mean = arrayfun(@(s) mean(s.summary.IDR_w), P)';
        dP95_worst_max = arrayfun(@(s) max(s.summary.dP95_worst), P)';
        T_end_worst_max = arrayfun(@(s) max(s.summary.T_end_worst), P)';
        idx_tbl = table(pol, ord, cd, qc_rate, PFA_w_mean, IDR_w_mean, dP95_worst_max, T_end_worst_max, ...
            'VariableNames',{'policy','order','cooldown_s','qc_rate','PFA_w_mean','IDR_w_mean','dP95_worst_max','T_end_worst_max'});
        writetable(idx_tbl, fullfile(outdir,'policy_index.csv'));
        for i = 1:numel(P)
            fname = sprintf('policy_%s_%s_cd%s.csv', ...
                Utils.sanitize_name(P(i).policy), Utils.sanitize_name(P(i).order), ...
                Utils.sanitize_name(num2str(P(i).cooldown_s)));
            writetable(P(i).summary, fullfile(outdir,fname));
        end
    catch
    end
end

end % end main function export_results

%% --- Helpers ---
function Cth = compute_Cth_effective(params)
    nStories = size(params.M,1) - 1;
    Rvec = params.toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
    mask = params.story_mask(:);  if numel(mask)==1, mask = mask*ones(nStories,1); end
    ndps = params.n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi = (mask .* ndps);
    V_oil_per = params.resFactor * (params.Ap * (2*params.Lgap));
    m_oil_tot = sum(multi) * (params.rho * V_oil_per);
    m_steel_tot = params.steel_to_oil_mass_ratio * m_oil_tot;
    Cth = max(m_oil_tot*params.cp_oil + m_steel_tot*params.cp_steel, eps);
end
