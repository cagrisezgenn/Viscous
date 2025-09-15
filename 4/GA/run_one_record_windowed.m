%% Kayıt Analiz Fonksiyonu
function out = run_one_record_windowed(rec, params, opts, prev_ts)
%RUN_ONE_RECORD_WINDOWED Pencereli metriklerle tek yer hareketini analiz eder.
%   OUT = RUN_ONE_RECORD_WINDOWED(REC, PARAMS, OPTS, PREV_TS) ölçekli kayıt
%   REC için doğrusal olmayan sönümleyici modeli çalıştırır ve Arias yoğunluğu
%   tabanlı bir zaman penceresi içinde performans metriklerini hesaplar. PARAMS
%   yapısı, yapı ve sönümleyici özelliklerini (bkz. PARAMETRELER.M) içerir.
%   OPTS davranışı kontrol eder:
%       window        - MAKE_ARIAS_WINDOW fonksiyonuna iletilen ayarlar
%       thermal_reset - 'each', 'carry' veya 'cooldown'
%       cooldown_s    - 'cooldown' modu için soğuma süresi [s]
%   PREV_TS, kayıtlara arasında termal durumu taşımak için isteğe bağlı bir
%   zaman serisi yapısına sahiptir.
%
%   OUT yapısı; name, scale, SaT1, win, metr, ts, qc_pass, T_start,
%   T_end, mu_end, clamp_hits, PFA, IDR, dP95,
%   Qcap95, cav_pct, t5, t95 ve coverage alanlarını içerir. TS
%   yalnızca COMPUTE_METRICS_WINDOWED için gerekli zaman serilerini barındırır.
%
%   Bu fonksiyon MCK_WITH_DAMPER ve COMPUTE_METRICS_WINDOWED
%   fonksiyonlarını kullanır.

% varsayılan argümanlar
if nargin < 4, prev_ts = []; end
if nargin < 3 || isempty(opts), opts = struct(); end

if isfield(opts,'thermal_reset') && strcmpi(opts.thermal_reset,'cooldown')
    if ~isfield(opts,'cooldown_s') || isempty(opts.cooldown_s) || isnan(opts.cooldown_s)
        opts.cooldown_s = 60;
    end
    opts.cooldown_s = max(opts.cooldown_s,0);
end

% QC eşikleri (eksik alanlar varsayılanlarla doldurulur)
if ~isfield(opts,'thr'), opts.thr = struct(); end
thr = Utils.default_qc_thresholds(opts.thr);

assert(isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), ...
    'run_one_record_windowed: params.thermal.hA_W_perK eksik');

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
catch
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
C_th = Utils.compute_Cth_effective(params);

switch mode
    case 'each'
        Tinit = params.T0_C;
    case 'carry'
        if ~isempty(prev_ts) && isfield(prev_ts,'T_oil')
            Tinit = prev_ts.T_oil(end);
        else
            Tinit = params.T0_C;
        end
    case 'cooldown'
        if ~isempty(prev_ts) && isfield(prev_ts,'T_oil')
            Tprev = prev_ts.T_oil(end);
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

%% Sönümleyicisiz Çözüm
% Lineer MCK sistemi sönümleyici olmadan çözülür.
[x0, a_rel0] = Utils.lin_MCK(rec.t, rec.ag, params.M, params.C0, params.K); %#ok<NASGU>

%% Damperli Çözüm
% Damperli çözüm mck_with_damper fonksiyonu üzerinden yürütülür.

mu_ref_eff   = params.mu_ref;
c_lam0_eff   = params.c_lam0;
try
    if isfield(params,'cfg') && isfield(params.cfg,'on') && isfield(params.cfg.on,'mu_floor') && params.cfg.on.mu_floor
        mu_min_phys = NaN;
        if isfield(params.cfg,'num') && isfield(params.cfg.num,'mu_min_phys') && isfinite(params.cfg.num.mu_min_phys)
            mu_min_phys = params.cfg.num.mu_min_phys;
        end
        if isfinite(mu_min_phys) && (mu_ref_eff < mu_min_phys)
            mu_ref_eff = mu_min_phys;
            scale_mu = mu_ref_eff / max(params.mu_ref, eps);
            c_lam0_eff = params.c_lam0 * scale_mu;
        end
    end
catch
    % sessiz geç
end
    [x,a_rel,ts] = mck_with_damper(rec.t, rec.ag, params.M, params.C0, params.K, ...
    params.k_sd, c_lam0_eff, params.Lori, params.orf, params.rho, params.Ap, ...
    params.Ao, params.Qcap_big, mu_ref_eff, params.thermal, Tinit, ...
    params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, params.Lgap, ...
    params.cp_oil, params.cp_steel, params.steel_to_oil_mass_ratio, ...
    params.story_mask, params.n_dampers_per_story, params.resFactor, params.cfg);

metr = compute_metrics_windowed(rec.t, x, a_rel, rec.ag, ts, params.story_height, win, params);

% Nominal koşu için son termal ve viskozite durumları kaydedilir
T_start = Tinit;
if isfield(ts,'T_oil')
    T_end = ts.T_oil(end);
    Tmax = params.T0_C + params.thermal.dT_max;
    clamp_hits = sum(diff(ts.T_oil >= Tmax) > 0);
else
    T_end = NaN;
    clamp_hits = NaN;
end
if isfield(ts,'mu')
    mu_end = ts.mu(end);
else
    mu_end = NaN;
end

qc_pass = (metr.cav_pct <= thr.cav_pct_max) && ...
          (metr.dP95 <= thr.dP95_max) && ...
          (metr.Qcap95 <= thr.Qcap95_max) && ...
          (T_end <= thr.T_end_max) && ...
          (mu_end >= thr.mu_end_min);

%% Çıktıların Derlenmesi
out = struct('name', rec.name, ...
    'scale', Utils.getfield_default(rec,'scale',1), ...
    'SaT1', Utils.getfield_default(rec,'IM',NaN), ...
    'win', win, 'metr', metr, 'ts', ts, 'qc_pass', qc_pass, ...
    'T_start', T_start, 'T_end', T_end, 'mu_end', mu_end, ...
    'clamp_hits', clamp_hits, ...
    'PFA', metr.PFA, 'IDR', metr.IDR, 'dP95', metr.dP95, 'Qcap95', metr.Qcap95, ...
    'cav_pct', metr.cav_pct, 't5', win.t5, 't95', win.t95, 'coverage', win.coverage);
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
catch
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




