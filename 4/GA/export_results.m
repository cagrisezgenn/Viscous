function export_results(outdir, scaled, params, opts, summary, all_out, varargin)
%EXPORT_RESULTS Sonuçları out/<timestamp>/ dizinine kaydeder.
%   EXPORT_RESULTS(OUTDIR, SCALED, PARAMS, OPTS, SUMMARY, ALL_OUT)
%   çeşitli günlükleri ve kayıt bazlı ayrıntıları yazar.
%   Ek yapılar (örneğin politika sonuçları) VARARGIN ile aktarılabilir.

if nargin < 1 || isempty(outdir), return; end
if ~exist(outdir,'dir'), mkdir(outdir); end

% Türetilmiş damper sabitlerini güncelle
params = Utils.recompute_damper_params(params);

%% Meta (IM alanlarını varsayılanlarla doldurma kaldırıldı)
    % TRIM bilgisi opsiyonel
    TRIM_names = Utils.getfield_default(opts,'TRIM_names',{});
    if isempty(TRIM_names), TRIM_names = {}; end

    % Varsa türetilmiş parametreler
    params_derived = struct();
        % Thermal mass derived parameter uses shared utility
        params_derived.C_th = Utils.compute_Cth_effective(params);

    % qc özet bilgisi
        pass = sum(summary.table.qc_pass);
        fail = height(summary.table) - pass;
        kshow = min(3, height(summary.table));
        [~,ix] = maxk(summary.table.PFA, kshow);
        worst3 = summary.table.name(ix);
        qc_meta = struct('pass',pass,'fail',fail,'worst3',{worst3});

    % snapshot.mat kaydı (P varsa dahil et) — IM_mode/band_fac/s_bounds yazılmaz
    snap = struct('params',params,'opts',opts,'scaled',scaled, ...
                  'TRIM_names',{TRIM_names},'params_derived',params_derived,'qc_meta',qc_meta);
    if ~isempty(varargin)
        snap.P = varargin{1};
    end
    save(fullfile(outdir,'snapshot.mat'), '-struct', 'snap', '-v7.3');

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
% ölçeklenmiş kayıt özetini CSV olarak yaz
writetable(tbl, fullfile(outdir,'scaled_index.csv'));

%% Özetler
% Genel özet tablosunu kaydet (zeta1_hot ve z2_over_z1_hot dahil)
writetable(summary.table, fullfile(outdir,'summary.csv'));

% Detaylı özet verisini .mat olarak sakla
full_summary = struct('summary',summary,'all_out',{all_out});
save(fullfile(outdir,'summary_full.mat'), '-struct', 'full_summary', '-v7.3');

% Kalite kontrol sonuçlarını yaz
pass = sum(summary.table.qc_pass);
fail = height(summary.table) - pass;
[~,idx] = maxk(summary.table.PFA, min(3,height(summary.table)));
worst = summary.table.name(idx);
qc_tbl = table(pass, fail, worst, 'VariableNames',{'pass','fail','worst'});
writetable(qc_tbl, fullfile(outdir,'qc_summary.csv'));
%% Kayıt Bazlı
% Her kayıt için ayrı klasör ve detay dosyaları oluştur
for k = 1:numel(all_out)
    out = all_out{k};
        recdir = fullfile(outdir, ['rec-' Utils.sanitize_name(out.name)]);
        if ~exist(recdir,'dir'), mkdir(recdir); end
        % pencere bilgisi (window.json)
        w = out.win;
        win_struct = struct('t5',w.t5,'t95',w.t95,'pad',Utils.getfield_default(w,'pad',0), ...
                            'coverage',w.coverage);
        Utils.writejson(win_struct, fullfile(recdir,'window.json'));
        % pencereye ait metrikler (metrics_win.csv)
        m = out.metr;
        mu_mode = Utils.getfield_default(out,'mu_mode','nominal');   % sadece etiket için
        mu_used = Utils.getfield_default(out,'mu_used',1.0);
        mstruct = struct('PFA_top',m.PFA_top,'IDR_max',m.IDR_max, ...
                         'dP95',m.dP_orf_q95,'Qcap95',m.Qcap_ratio_q95, ...
                         'cav_pct',m.cav_pct,'T_end',out.T_end,'mu_end',out.mu_end, ...
                         'P_mech',Utils.getfield_default(m,'P_mech',NaN), ...
                         'Re_max',Utils.getfield_default(m,'Re_max',NaN), ...
                         'zeta1_hot',Utils.getfield_default(m,'zeta1_hot',NaN), ...
                         'z2_over_z1_hot',Utils.getfield_default(m,'z2_over_z1_hot',NaN));
        % OUT yapısında PF bilgisi varsa ekle
        if isfield(out,'PF_t_on'),       mstruct.PF_t_on = out.PF_t_on; end
        if isfield(out,'PF_tau'),        mstruct.PF_tau = out.PF_tau; end
        if isfield(out,'PF_gain'),       mstruct.PF_gain = out.PF_gain; end
        if isfield(out,'PF_mode'),       mstruct.PF_mode = out.PF_mode; end
        if isfield(out,'PF_auto_t_on'),  mstruct.PF_auto_t_on = out.PF_auto_t_on; end
        mstruct.mu_mode = mu_mode; mstruct.mu_used = mu_used;
        % Yeni parametreleri isteğe bağlı olarak metriklere ekle
        param_fields = {'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref'};
        for ii = 1:numel(param_fields)
            fn = param_fields{ii};
            if isfield(out, fn)
                mstruct.(fn) = out.(fn);
            end
        end
        if isfield(m,'E_orifice_win'), mstruct.E_orf_win = m.E_orifice_win; end
        if isfield(m,'E_struct_win'), mstruct.E_struct_win = m.E_struct_win; end
        metr_tbl = struct2table(mstruct);
        writetable(metr_tbl, fullfile(recdir,'metrics_win.csv'));
        % seyreltilmiş zaman serileri (ts_ds.mat)
        if isfield(out,'ts') && ~isempty(out.ts)
                ds = 5;
                if isfield(opts,'export') && isfield(opts.export,'ds')
                    ds = opts.export.ds;
                end
                ts_ds = Utils.downsample_ts(out.ts, ds);
                tmp_struct = struct('ts_ds',ts_ds);
                save(fullfile(recdir,'ts_ds.mat'), '-struct', 'tmp_struct', '-v7.3');
        end
        % isteğe bağlı grafikler
        if isfield(opts,'export') && isfield(opts.export,'plots') && opts.export.plots
                figdir = fullfile(recdir,'figures');
                if ~exist(figdir,'dir'), mkdir(figdir); end
                f = figure('Visible','off');
                bar(out.metr.IDR_max); title('IDR_{max}');
                saveas(f, fullfile(figdir,'IDR_max.png')); close(f);
                f = figure('Visible','off');
                bar(out.metr.PFA_top); title('PFA_{top}');
                saveas(f, fullfile(figdir,'PFA_top.png')); close(f);
        end
end

%%% Politika Sonuçları
% Politika sonuçlarını özetleyip kaydet
if ~isempty(varargin)
    P = varargin{1};
    pol = {P.policy}';
    ord = {P.order}';
    cd  = [P.cooldown_s]';
    qc_rate = arrayfun(@(s) s.qc.pass_fraction, P)';
    PFA_mean = arrayfun(@(s) mean(s.summary.PFA), P)';
    IDR_mean = arrayfun(@(s) mean(s.summary.IDR), P)';
    dP95_max = arrayfun(@(s) max(s.summary.dP95), P)';
    T_end_max = arrayfun(@(s) max(s.summary.T_end), P)';
    idx_tbl = table(pol, ord, cd, qc_rate, PFA_mean, IDR_mean, dP95_max, T_end_max, ...
        'VariableNames',{'policy','order','cooldown_s','qc_rate','PFA_mean','IDR_mean','dP95_max','T_end_max'});
    writetable(idx_tbl, fullfile(outdir,'policy_index.csv'));
    for i = 1:numel(P)
        fname = sprintf('policy_%s_%s_cd%s.csv', ...
            Utils.sanitize_name(P(i).policy), Utils.sanitize_name(P(i).order), ...
            Utils.sanitize_name(num2str(P(i).cooldown_s)));
        writetable(P(i).summary, fullfile(outdir,fname));
    end
end

end % export_results fonksiyon sonu

