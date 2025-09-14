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
try
    % TRIM bilgisi opsiyonel
    TRIM_names = Utils.getfield_default(opts,'TRIM_names',{});
    if isempty(TRIM_names), TRIM_names = {}; end

    % Varsa türetilmiş parametreler
    params_derived = struct();
    try
        params_derived.C_th = compute_Cth_effective(params);
    catch
    end

    % qc özet bilgisi
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

    % snapshot.mat kaydı (P varsa dahil et) — IM_mode/band_fac/s_bounds yazılmaz
    snap = struct('params',params,'opts',opts,'scaled',scaled, ...
                  'TRIM_names',{TRIM_names},'params_derived',params_derived,'qc_meta',qc_meta);
    if ~isempty(varargin)
        snap.P = varargin{1};
    end
    safe_write(snap, fullfile(outdir,'snapshot.mat'), @(d,f) save(f,'-struct','d','-v7.3'));
catch
end

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
safe_write(tbl, fullfile(outdir,'scaled_index.csv'), @writetable);

%% Özetler
% Genel özet tablosunu kaydet (zeta1_hot ve z2_over_z1_hot dahil)
safe_write(summary.table, fullfile(outdir,'summary.csv'), @writetable);

% Detaylı özet verisini .mat olarak sakla
full_summary = struct('summary',summary,'all_out',{all_out});
safe_write(full_summary, fullfile(outdir,'summary_full.mat'), @(d,f) save(f,'-struct','d','-v7.3'));

% Kalite kontrol sonuçlarını yaz
pass = sum(summary.table.qc_all_mu);
fail = height(summary.table) - pass;
[~,idx] = maxk(summary.table.PFA_worst, min(3,height(summary.table)));
worst = summary.table.name(idx);
qc_tbl = table(pass, fail, worst, 'VariableNames',{'pass','fail','worst'});
safe_write(qc_tbl, fullfile(outdir,'qc_summary.csv'), @writetable);
%% Kayıt Bazlı
% Her kayıt için ayrı klasör ve detay dosyaları oluştur
for k = 1:numel(all_out)
    out = all_out{k};
    try
        recdir = fullfile(outdir, ['rec-' Utils.sanitize_name(out.name)]);
        if ~exist(recdir,'dir'), mkdir(recdir); end
        % pencere bilgisi (window.json)
        w = out.win;
        win_struct = struct('t5',w.t5,'t95',w.t95,'pad',Utils.getfield_default(w,'pad',0), ...
                            'coverage',w.coverage);
        safe_write(win_struct, fullfile(recdir,'window.json'), @Utils.writejson);
        % mu sonuçları (mu_results.csv)
        if isfield(out,'mu_results') && ~isempty(out.mu_results)
            mu_tbl = struct2table(arrayfun(@(s) s.metr, out.mu_results));
            mu_tbl.mu = [out.mu_results.mu_factor]';
            mu_tbl = movevars(mu_tbl,'mu','before',1);
            safe_write(mu_tbl, fullfile(recdir,'mu_results.csv'), @writetable);
        end
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
        safe_write(metr_tbl, fullfile(recdir,'metrics_win.csv'), @writetable);
        % seyreltilmiş zaman serileri (ts_ds.mat)
        if isfield(out,'ts') && ~isempty(out.ts)
            try
                ds = 5;
                if isfield(opts,'export') && isfield(opts.export,'ds')
                    ds = opts.export.ds;
                end
                ts_ds = Utils.downsample_ts(out.ts, ds);
                safe_write(struct('ts_ds',ts_ds), fullfile(recdir,'ts_ds.mat'), @(d,f) save(f,'-struct','d','-v7.3'));
            catch
            end
        end
        % isteğe bağlı grafikler
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

%%% Politika Sonuçları
% Politika sonuçlarını özetleyip kaydet
if ~isempty(varargin)
    P = varargin{1};
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
    safe_write(idx_tbl, fullfile(outdir,'policy_index.csv'), @writetable);
    for i = 1:numel(P)
        fname = sprintf('policy_%s_%s_cd%s.csv', ...
            Utils.sanitize_name(P(i).policy), Utils.sanitize_name(P(i).order), ...
            Utils.sanitize_name(num2str(P(i).cooldown_s)));
        safe_write(P(i).summary, fullfile(outdir,fname), @writetable);
    end
end

end % export_results fonksiyon sonu

%%% --- Yardımcı Fonksiyonlar ---
function Cth = compute_Cth_effective(params)
    % Sıcaklık kapasitesinin etkin değerini hesaplar
    nStories = size(params.M,1) - 1;
        mask = params.story_mask(:);  if numel(mask)==1, mask = mask*ones(nStories,1); end
    ndps = params.n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi = (mask .* ndps);
    V_oil_per = params.resFactor * (params.Ap * (2*params.Lgap));
    m_oil_tot = sum(multi) * (params.rho * V_oil_per);
    m_steel_tot = params.steel_to_oil_mass_ratio * m_oil_tot;
    Cth = max(m_oil_tot*params.cp_oil + m_steel_tot*params.cp_steel, eps);
end

function safe_write(obj, filepath, writeFcn)
% Verilen yazma fonksiyonunu güvenle çağırır, hataları uyarı olarak bildirir
    try
        writeFcn(obj, filepath);
    catch ME
        warning('Yazma hatası (%s): %s', filepath, ME.message);
    end
end

