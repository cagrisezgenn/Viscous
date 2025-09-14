%% Kayıt Analiz Fonksiyonu
function out = run_one_record_windowed(rec, params, opts, prev_diag)
%RUN_ONE_RECORD_WINDOWED Pencereli metriklerle tek yer hareketini analiz eder.
%   OUT = RUN_ONE_RECORD_WINDOWED(REC, PARAMS, OPTS, PREV_DIAG) ölçekli kayıt
%   REC için doğrusal olmayan sönümleyici modelini çalıştırır ve Arias yoğunluğu
%   tabanlı bir zaman penceresi içinde performans metriklerini hesaplar. PARAMS
%   yapısı, yapı ve sönümleyici özelliklerini (bkz. PARAMETRELER.M) içerir.
%   OPTS davranışı kontrol eder:
%       window        - MAKE_ARIAS_WINDOW fonksiyonuna iletilen ayarlar
%       thermal_reset - 'each', 'carry' veya 'cooldown'
%       cooldown_s    - 'cooldown' modu için soğuma süresi [s]
%   PREV_DIAG, kayıtlara arasında termal durumu taşımak için isteğe bağlı bir
%   diagnostiğe sahiptir.
%
%   OUT yapısı; name, scale, SaT1, win, metr, diag, mu_results, weighted,
%   worst, ts, qc_all_mu, T_start, T_end, mu_end, clamp_hits, PFA_top,
%   IDR_max, dP_orf_q95, Qcap_ratio_q95, cav_pct, t5, t95 ve coverage
%   alanlarını içerir.
%
%   Bu fonksiyon, dosya sonunda yer alan MCK_WITH_DAMPER_TS yardımcı
%   rutini ve COMPUTE_METRICS_WINDOWED fonksiyonunu kullanır.

% varsayılan argümanlar
if nargin < 4, prev_diag = []; end
if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts,'mu_factors'), opts.mu_factors = 1; end
if ~isfield(opts,'mu_weights'), opts.mu_weights = 1; end

% Türetilmiş damper sabitlerini güncelle
params = Utils.recompute_damper_params(params);

if isfield(opts,'thermal_reset') && strcmpi(opts.thermal_reset,'cooldown')
    if ~isfield(opts,'cooldown_s') || isempty(opts.cooldown_s) || isnan(opts.cooldown_s)
        opts.cooldown_s = 60;
    end
    opts.cooldown_s = max(opts.cooldown_s,0);
end

% QC eşikleri (eksik alanlar varsayılanlarla doldurulur)
if ~isfield(opts,'thr'), opts.thr = struct(); end
thr = Utils.default_qc_thresholds(opts.thr);
opts.thr = thr;

assert(isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), ...
    'run_one_record_windowed: params.thermal.hA_W_perK eksik');

assert(numel(opts.mu_factors)==numel(opts.mu_weights), ...
    'mu_factors and mu_weights must have same length.');
mu_weights = opts.mu_weights(:);
wsum = sum(mu_weights);
assert(wsum>0,'mu_weights sum must be > 0.');
mu_weights = mu_weights/wsum;
mu_factors = opts.mu_factors(:)';

%% Pencere Hazırlığı
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

% PF auto_t_on, Arias t5'e göre belirlenir (çözücü çağrısından önce)
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
catch ME
    warning('PF auto_t_on başarısız: %s', ME.message);
    % parametreler değiştirilmeden bırakıldı
end

% ham ve ölçekli kayıt karşılaştırması kaldırıldı (üst seviyede kullanılmıyor)

%% Termal Sıfırlama
Tinit = params.T0_C;
if isfield(opts,'thermal_reset')
    mode = opts.thermal_reset;
else
    mode = 'each';
end

% soğuma seçeneği için termal kapasite hesaplanır
 nStories = size(params.M,1) - 1;
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

%% Sönümleyicisiz Çözüm
% Lineer MCK sistemi sönümleyici olmadan çözülür.
[x0, a_rel0] = Utils.lin_MCK(rec.t, rec.ag, params.M, params.C0, params.K); %#ok<NASGU>

%% Damperli Çözüm
% Damperli çözüm doğrudan mck_with_damper_ts fonksiyonu üzerinden yürütülür.

nMu = numel(mu_factors);
mu_results = struct('mu_factor',cell(1,nMu));

for i = 1:nMu
    f = mu_factors(i);
    mu_ref_eff   = params.mu_ref  * f;
    c_lam0_eff   = params.c_lam0  * f;
    % Adım 0: Yüksek sıcaklıkta μ düşüşüne karşı taban (opsiyonel)
    try
        if isfield(params,'cfg') && isfield(params.cfg,'on') && isfield(params.cfg.on,'mu_floor') && params.cfg.on.mu_floor
            mu_min_phys = NaN;
            if isfield(params.cfg,'num') && isfield(params.cfg.num,'mu_min_phys') && isfinite(params.cfg.num.mu_min_phys)
                mu_min_phys = params.cfg.num.mu_min_phys;
            end
            if isfinite(mu_min_phys) && (mu_ref_eff < mu_min_phys)
                mu_ref_eff = mu_min_phys;
                % c_lam0 ~ μ ile orantılı olduğundan, tutarlılık için yeniden ölçekle
                scale_mu = mu_ref_eff / max(params.mu_ref, eps);
                c_lam0_eff = params.c_lam0 * scale_mu;
            end
        end
    catch
        % sessiz geç
    end
    [x,a_rel,ts,diag] = mck_with_damper_ts(rec.t, rec.ag, params.M, params.C0, params.K, ...
        params.k_sd, c_lam0_eff, params.Lori, opts.use_orifice, params.orf, params.rho, params.Ap, ...
        params.A_o, params.Qcap_big, mu_ref_eff, opts.use_thermal, params.thermal, ...
        params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, params.Lgap, ...
        params.cp_oil, params.cp_steel, params.steel_to_oil_mass_ratio, ...
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

% Nominal metrikler (f=1)
[~,nom_idx] = min(abs(mu_factors-1));
metr = mu_results(nom_idx).metr;
diag = mu_results(nom_idx).diag;

% Nominal koşu için son termal ve viskozite durumları kaydedilir
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

% Ağırlıklı ve en kötü durum özetleri (metrik bazında uygun min/maks)
fields = {'PFA_top','IDR_max','dP_orf_q95','dP_orf_q50','Q_q95','Q_q50','Qcap_ratio_q95', ...
          'cav_pct','T_oil_end','mu_end', 'x10_max_D','a10abs_max_D', ...
          'E_orifice_full','E_struct_full','E_ratio_full','PF_p95'};
weighted = struct();
worst = struct();
worst.which_mu = struct();
for kf = 1:numel(fields)
    fn = fields{kf};
    vals = arrayfun(@(s) s.metr.(fn), mu_results);
    % tüm metrikler için ağırlıklı ortalama
    weighted.(fn) = sum(mu_weights(:)'.*vals);
    % en kötü durum seçim kuralı
    switch fn
        case 'mu_end'
            [worst.(fn), idx] = min(vals); % daha küçük viskozite daha kötüdür
        otherwise
            [worst.(fn), idx] = max(vals);
    end
    worst.which_mu.(fn) = mu_results(idx).mu_factor;
end

%% Çıktıların Derlenmesi
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
% kolaylık amaçlı telemetri alanları
out.PFA_top = metr.PFA_top;
out.IDR_max = metr.IDR_max;
out.dP_orf_q95 = metr.dP_orf_q95;
out.Qcap_ratio_q95 = metr.Qcap_ratio_q95;
out.cav_pct = metr.cav_pct;
out.t5 = win.t5; out.t95 = win.t95; out.coverage = win.coverage;
% ayarlanmışsa PF rampasının başlangıcı
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
catch ME
    warning('PF telemetri yakalama başarısız: %s', ME.message);
end

% Yeni parametreleri isteğe bağlı olarak kaydet
param_fields = {'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref'};
for ii = 1:numel(param_fields)
    fn = param_fields{ii};
    if isfield(params, fn)
        out.(fn) = params.(fn);
    end
end
end

%% =====================================================================
%% mck_with_damper_ts Yerel Fonksiyonu
%% =====================================================================
function [x,a_rel,ts,diag] = mck_with_damper_ts(t,ag,M,C,K, k_sd,c_lam0,Lori, ...
    use_orifice, orf, rho, Ap, Ao, Qcap, mu_ref, use_thermal, thermal, ...
    T_ref_C, b_mu, c_lam_min, c_lam_cap, Lgap, cp_oil, cp_steel, ...
    steel_to_oil_mass_ratio, story_mask, n_dampers_per_story, ...
    resFactor, cfg)
%1) MCK_WITH_DAMPER çözümü
[x,a_rel,diag] = mck_with_damper(t,ag,M,C,K, k_sd,c_lam0,Lori, use_orifice, orf, ...
    rho, Ap, Ao, Qcap, mu_ref, use_thermal, thermal, T_ref_C, b_mu, ...
    c_lam_min, c_lam_cap, Lgap, cp_oil, cp_steel, steel_to_oil_mass_ratio, ...
    story_mask, n_dampers_per_story, resFactor, cfg);

%2) Öykü vektörlerinin hazırlanması
nStories = size(diag.drift,2);
mask = story_mask(:); if numel(mask)==1, mask = mask*ones(nStories,1); end
ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
multi = (mask .* ndps).';

%3) Zaman serisi yapılarının oluşturulması
ts = struct();
ts.t = t;
ts.drift = diag.drift;
ts.dvel = diag.dvel;
ts.story_force = diag.story_force;
ts.PF = diag.PF;
ts.Q = diag.Q;
ts.dP_orf = diag.dP_orf;
ts.Qcap_ratio = abs(diag.Q) ./ Qcap;
ts.cav_mask = diag.dP_orf < 0;

%4) Güç bileşenlerinin hesabı
P_visc_per = diag.c_lam .* (diag.dvel.^2);
ts.P_visc = sum(P_visc_per .* multi, 2);

if isfield(diag,'P_orf_per') && ~isempty(diag.P_orf_per)
    P_orf_per = diag.P_orf_per;
else
    P_orf_per = abs(diag.dP_orf .* diag.Q);
end
ts.P_orf = sum(P_orf_per .* multi, 2);

if isfield(diag,'P_sum')
    ts.P_sum = diag.P_sum;
elseif isfield(ts,'P_orf') && isfield(ts,'P_visc')
    ts.P_sum = ts.P_orf + ts.P_visc;
else
    ts.P_sum = [];
end

%5) Enerji birikimlerinin hesaplanması
P_struct = sum(diag.story_force .* diag.dvel, 2);
ts.E_orf = cumtrapz(t, ts.P_orf);
ts.E_struct = cumtrapz(t, P_struct);

% DIAG yapısı değiştirilmeden geri döndürülür
end




