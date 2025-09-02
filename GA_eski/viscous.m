%% =====================================================================
%  10-KatlÄ± Ã‡erÃ§eve â€” ODE-only viskoz damper modeli (CLEAN / AHENK + Switch)
%  (SÄ±kÄ±ÅŸabilirlik + Hat ataletÄ± + Cd(Re) orifis + kavitasyon + 2-dÃ¼ÄŸÃ¼mlÃ¼ Ä±sÄ±)
%  Laminer kayÄ±p hidrolikte: Î”p_lam(T) = R_lam(T)*Q
%  Kavitasyon: p2_eff = max(p2, cav_sf * p_vap(T)) hem akÄ±ÅŸta hem kuvvette
%  Solver: ode23tb (+ gÃ¼venli deval + ode15s fallback)
%  NOT: Antoine katsayÄ±larÄ±nÄ± yaÄŸÄ±nÄ±za gÃ¶re kalibre edin.
% =====================================================================

clear; clc; close all; clear ga;
[paths, env] = init_environment();

%% -------------------- KullanÄ±cÄ± anahtarlarÄ± ---------------------------

% (1) KayÄ±t/yÃ¶n seÃ§imi + gÃ¶rselleÅŸtirme
vis.make_plots          = true;
vis.sel_rec             = 1;         % 1..R
vis.sel_dir             = 'x';       % 'X' veya 'Y'
vis.dual_run_if_needed  = true;      % plot_source â‰  sim_source ise ikinci koÅŸu yap

% (2) Ã–rnekleme / yeniden Ã¶rnekleme anahtarlarÄ±
prep.target_dt      = 0.005;     % hedef dt (s)
prep.resample_mode  = 'auto';     % 'auto'|'off'|'force'
prep.regrid_method  = 'pchip';   % 'pchip'|'linear'
prep.tol_rel        = 1e-6;      % dt eÅŸitlik toleransÄ± (gÃ¶reli)

% (3) Åiddet eÅŸitleme / PSA / CMS anahtarlarÄ±
pp.on.intensity      = true;       % band-ortalama Sa(+CMS) normalizasyonu (yalnÄ±z GA iÃ§in)
pp.on.CMS            = false;      % cms_target.mat varsa ve kullanmak istersen true
pp.gammaCMS          = 0.50;       % hibrit aÄŸÄ±rlÄ±k (0â†’yalnÄ±z band, 1â†’yalnÄ±z CMS)

pp.PSA.zeta          = 0.05;       % SDOF sÃ¶nÃ¼m oranÄ±
pp.PSA.band_fac      = [0.8 1.2];  % T1 bandÄ± (T1Â±%20)
pp.PSA.Np_band       = 21;         % band iÃ§i periyot sayÄ±sÄ±
pp.PSA.downsample_enabled = true;  % SA hesabÄ± iÃ§in downsample kullan (false â†’ ham dt)
pp.PSA.downsample_dt      = 0.02;  % [s] downsample hedef dt; yalnÄ±z downsample_enabled=true iken
pp.PSA.use_parfor        = true;  % Parallel Toolbox varsa denersin

% (3b) SimÃ¼lasyon paralelleÅŸtirme anahtarÄ±
sim.use_parfor           = false;  % Âµ-senaryolarÄ±nÄ± parfor ile Ã§alÄ±ÅŸtÄ±r

% (4) Arias penceresi & kuyruk
pp.on.arias          = true;
pp.tail_sec          = 10;

% (5) SimÃ¼lasyon ve grafik iÃ§in veri kaynaÄŸÄ± seÃ§imi
%     'raw'    : ham ivmeler (genlik korunur)
%     'scaled' : ÅŸiddet eÅŸitlenmiÅŸ ivmeler (GA/normalize amaÃ§lÄ±)
sel.sim_source  = 'scaled';
sel.plot_source = 'scaled';

%% -------------------- AmaÃ§ fonksiyonu anahtarlarÄ± ---------------------

% (A) Hedef katlar ve metrik
obj.idx_disp_story   = 10;        % d_10 â†’ tepe yer deÄŸiÅŸtirme
obj.idx_acc_story    = 10;         % a_3  â†’ ivme metriÄŸi
obj.acc_metric       = 'rms+p95';     % 'rms' | 'energy' | 'rms+p95' (hibrit)
obj.p95_penalty_w    = 0.40;      % hibritte pik cezasÄ± aÄŸÄ±rlÄ±ÄŸÄ±

%% -------------------- KÄ±sÄ±t anahtarlarÄ± -------------------------------

% KÄ±sÄ±tlarÄ± aÃ§/kapa ve eÅŸik/ceza ayarlarÄ±
cons.on.spring_tau     = false;    % K1: yay kesme gerilmesi (Ï„_max â‰¤ Ï„_allow)
cons.on.spring_slender = false;    % K2: L_free/D_m â‰¤ Î»_max
cons.on.stroke         = false;    % K3: max|drift| â‰¤ 0.9*L_gap
cons.on.force_cap      = false;    % K4: max|F_story| â‰¤ F_cap
cons.on.dp_quant       = false;    % K5: qâ‰ˆ0.99 Î”p_orf â‰¤ dP_cap
cons.on.thermal_dT     = false;    % K6: Î”T_est â‰¤ Î”T_cap
cons.on.cav_frac       = false;    % K7: kavitasyon payÄ± sÄ±nÄ±rÄ±
cons.on.qsat_margin    = false;    % K8: Q 95p â‰¤ margin*Qcap_big
cons.on.fail_bigM      = false;    % K9: SimÃ¼lasyon baÅŸarÄ±sÄ±zsa bÃ¼yÃ¼k ceza

% AÅŸama-Ã¶zel (set 1/2/3) kÄ±sÄ±t otomatik aÃ§ma kapama
% true: aÅŸaÄŸÄ±daki set-Ã¶zgÃ¼ bloklar kÄ±sÄ±tlarÄ± otomatik olarak aÃ§abilir (mevcut davranÄ±ÅŸ)
% false: kullanÄ±cÄ±daki Ã¼stteki cons.on.* ayarlarÄ± her aÅŸamada birebir korunur
cons.allow_stage_overrides = true;

% EÅŸikler / limitler
cons.spring.tau_allow   = 300e6;  % [Pa] yay Ã§eliÄŸi iÃ§in tipik, gerekirse gÃ¼ncelle
cons.spring.lambda_max  = 12.0;   % boy/Ã§ap sÄ±nÄ±rÄ±
% true -> serbest yay boyu sabit deÄŸerden alÄ±nÄ±r; false -> boÅŸluk oranÄ±na gÃ¶re hesaplanÄ±r
cons.spring.use_fixed_length = false;   % varsayÄ±lan: false
cons.spring.L_free_fixed     = NaN;     % [m] use_fixed_length=true ise kullanÄ±lan serbest boy
cons.spring.L_free_auto_fac = 2.2; % 'auto' modda L_free â‰ˆ fac * L_gap

cons.stroke.util_factor = 0.90;   % izinli strok = 0.90*L_gap

cons.force.F_cap        = 2e6;    % [N] cihaz kuvvet sÄ±nÄ±rÄ± (varsayÄ±lan: 2e6); kÄ±sÄ±tÄ± cons.on.force_cap ile aÃ§/kapa
cons.dp.q               = 0.99;   % Î”p_orf zaman-iÃ§i quantile
cons.dp.agg             = 'max';  % kayÄ±tlar arasÄ±: 'max' | 'cvar'
cons.alpha_CVaR_cons    = 0.20;   % yalnÄ±z 'cvar' seÃ§ilirse kullanÄ±lÄ±r

cons.thermal.cap_C      = 30.0;   % [Â°C] yaÄŸ Ä±sÄ±nma sÄ±nÄ±rÄ±
cons.hyd.cav_frac_cap   = 0.05;   % kavitasyon zaman oranÄ± sÄ±nÄ±rÄ± (95p)
cons.hyd.Q_margin       = 0.90;   % Q 95p â‰¤ margin*Qcap_big

% Ceza ayarlarÄ±
cons.pen.power  = 2;              % hinge^power
cons.pen.bigM   = 1e6;            % simÃ¼lasyon baÅŸarÄ±sÄ±zlÄ±ÄŸÄ±nda eklenecek ceza
cons.pen.lambda = struct( ...     % kÄ±sÄ±t baÅŸÄ±na aÄŸÄ±rlÄ±klar
    'spring_tau',     0, ...%0.3
    'spring_slender', 0, ...%0.2
    'stroke',         0, ...%1.2
    'force_cap',      0, ...%0
    'dp_quant',       0, ...%0.5
    'thermal_dT',     0, ...%0.2
    'cav_frac',       0, ...%1.5
    'qsat_margin',    0, ...%0.3
    'fail_bigM',      1 );%

% KÄ±sÄ±t simÃ¼lasyon kaynaÄŸÄ± ve Î¼ seti (amaÃ§la tutarlÄ± kalsÄ±n)
cons.src_for_constraints = 'raw';           % kÄ±sÄ±tlarÄ± 'raw' ile deÄŸerlendir
if isfield(obj,'mu_scenarios')
    cons.mu_scenarios = obj.mu_scenarios;
else
    cons.mu_scenarios = [0.75 1.00 1.25];  % varsayÄ±lan; aÅŸaÄŸÄ±da obj tanÄ±mlanÄ±nca senkronlayacaÄŸÄ±z
end

% (B) Arias penceresi iÃ§i Ã¶lÃ§Ã¼m
obj.use_arias_window = true;      % true â†’ [t5,t95] aralÄ±ÄŸÄ±nda metrik
obj.window_source    = 'same';    % 'same' â†’ kullanÄ±lan seri Ã¼zerinden hesapla

% (C) YÃ¶n zarfÄ±
obj.dir_mode         = 'Xonly'; % 'envelope' (=max(X,Y)) | 'Xonly' | 'Yonly'

% (D) Î¼-robustluk
obj.mu_scenarios = [0.9 1.0 1.1];
obj.mu_aggregate = 'weighted';
obj.mu_weights   = [0.25 0.50 0.25];
% KÄ±sÄ±t Î¼-senaryolarÄ±nÄ± amaÃ§la hizala
cons.mu_scenarios = obj.mu_scenarios;

% (E) Risk agregasyonu (kayÄ±tlar Ã¼zerinde)
obj.use_scaled_for_goal = true;   % amaÃ§ iÃ§in 'scaled' setiyle Ã§alÄ±ÅŸ (Ã¶nerilir)
obj.alpha_CVaR       = 0.20;      % tail payÄ±
obj.weights_da       = [0.5 0.5]; % [w_disp, w_acc] toplam 1.0

% (F) Referans tanÄ±mÄ±
% d_ref ve a_ref: aynÄ± kayÄ±t/yÃ¶n iÃ§in (damper yok, Î¼=1.0), sabit baz koÅŸusundan.
% NOT: Bu referanslar bir kez hesaplanÄ±r ve tasarÄ±mla deÄŸiÅŸmez.

%% -------------------- Model anahtarlarÄ± (blok-bazlÄ± aÃ§/kapa) ----------

% NOT: Bu blok dosyanÄ±n BAÅINDA olmalÄ±. AÅŸaÄŸÄ±daki varsayÄ±lanlarÄ±
% dilediÄŸin gibi deÄŸiÅŸtir; ayrÄ±ca ensure_cfg_defaults() eksikleri tamamlar.
cfg.use_orifice = true;
cfg.use_thermal = true;
cfg.on.CdRe            = true;
cfg.on.Rlam            = true;
cfg.on.Rkv             = true;
cfg.on.Qsat            = true;
cfg.on.cavitation      = true;
cfg.on.dP_cap          = true;
cfg.on.hyd_inertia     = true;
cfg.on.leak            = true;
cfg.on.pressure_ode    = true;
cfg.on.pressure_force  = true;
cfg.on.mu_floor        = true;


% BasÄ±nÃ§-kuvvet rampa/kazanÃ§ (PF)
cfg.PF.mode       = 'ramp';
cfg.PF.auto_t_on  = false;   % true -> t_on = t5 + 0.5 via set_pf_ton_if_auto
cfg.PF.t_on       = 2;      % auto_t_on=false ise manuel deÄŸer (s)
cfg.PF.tau        = 1.5;    % 4.0 â†’ 1.5  (Ã¶nerim)
cfg.PF.k         = 2;  % softness for PF ramp
cfg.PF.gain       = 0.6;  % PF gain within 0.4-0.8
cfg.on.pf_resistive_only = true;  % sadece rezistif PF bileÅŸeni kullan

% Eksik/yanlÄ±ÅŸ alanlar iÃ§in gÃ¼venli tamamlayÄ±cÄ± (guard; PF.auto_t_on dahil)
cfg = ensure_cfg_defaults(cfg);

%% -------------------- GA/Opt tasarÄ±m anahtarlarÄ± ----------------------

ga.enable      = false;    % GA tasarÄ±m vektÃ¶rÃ¼ uygula? (false â†’ baz set)
ga.design_set  = 1;        % 1|2|3 (aÅŸaÄŸÄ±daki set tanÄ±mlarÄ±)
ga.x           = [];       % Ã–rn: set-1 iÃ§in 9x1 vektÃ¶r, boÅŸsa uygulanmaz
% Ã–rnek: otomatik orta-nokta denemesi
% [lb,ub,~,~] = ga_get_bounds(ga.design_set); ga.x = 0.5*(lb+ub); ga.enable=true;

%% -------------------- GA Ã‡ALIÅTIRICI (opsiyonel) ----------------------

% --- DeÄŸiÅŸken sÄ±nÄ±rlarÄ± ve IntCon hazÄ±rlÄ±ÄŸÄ± ---
[lb,ub,int_idx,~] = ga_get_bounds(ga.design_set);
if isempty(int_idx)
    IntCon = [];           % tamsayÄ± deÄŸiÅŸken yok
else
    IntCon = int_idx(:)';  % tamsayÄ± indeks vektÃ¶rÃ¼
end

% GA anahtarlarÄ±
ga.opt.enable       = false;
ga.opt.seed         = 11;
ga.opt.popsize      = [48];
% Ã‡ok aÅŸamalÄ± akÄ±ÅŸ: Ã¶nce set-1, sonra set-3
ga.opt.multi_stage  = [2 2 1 1 5 5 3 3 4 4];
% AÅŸama baÅŸÄ±na nesil sayÄ±sÄ± (opsiyonel):
ga.opt.maxgen_stage = [20 10 20 10 16 8 16 8 12 6];  % her blok iÃ§in nesil
ga.opt.maxgen       = ga.opt.maxgen_stage(1);   % ilk aÅŸama iÃ§in baÅŸlangÄ±Ã§
ga.opt.use_parallel  = true;      % GA.UseParallel
ga.opt.lhs_refine    = true;      % 2-aÅŸamalÄ± LHS (coarseâ†’refine)
ga.opt.cache_enable  = true;      % aynÄ± x iÃ§in memoize
ga.opt.fail_early_k  = 10;         % erken Ã§Ä±kÄ±ÅŸ eÅŸiÄŸi (baÅŸarÄ±sÄ±z sim sayÄ±sÄ±)
ga.opt.save_log      = 'runLog_ga.mat';  % nÃ¼fus & en iyi Ã§Ã¶zÃ¼m logu
% daraltma parametreleri
ga.opt.keep_top     = [0.35 0.35 0.30 0.30 0.25 0.25 0.20 0.20 0.20 0.20];   % en iyi %20'yi tut
ga.opt.buffer_fac   = 0.10;   % p10â€“p90 etrafÄ±na %10 tampon

%% -------------------- Girdiler (7 kayÄ±t; sÃ¼tunlar: t, ax, (ops) ay) ---

Sall = load('acc_matrix.mat');
fn   = fieldnames(Sall);
fn   = fn(startsWith(fn,'acc_matrix'));
R    = numel(fn);
if R==0, error('acc_matrix.mat iÃ§inde acc_matrix* isimli dizi bulunamadÄ±.'); end

%% -------------------- YapÄ± (T1 iÃ§in gerekli) --------------------------

n  = 10;
m  = 360e3 * ones(n,1);
k  = 6.5e8 * ones(n,1);
c0  = 6.2e6 * ones(n,1);
[M,K,Cstr] = make_KCM(n,m,k,c0);
[~,D] = eig(K,M); w = sqrt(sort(diag(D),'ascend')); T1 = 2*pi/w(1);
% --- IDR iÃ§in kat yÃ¼kseklik(leri) ---
h_story_m = 3.0 * ones(n-1,1);   % tÃ¼m katlar 3.0 m ise
% h_story_m = [h1; h2; ...; h_{n-1}];   % kat kat farklÄ± ise (alternatif)


%% -------------------- PSA fonksiyon seÃ§imi ----------------------------

f_band = @sdof_PSA_band_avg_aug;   % tek ODE, Ã§ok periyot
f_vec  = @sdof_PSA_vec_aug_ode;

%% -------------------- KayÄ±tlarÄ± oku â†’ RAW (tekil & hedef dt kontrolÃ¼) -

t_rawX = cell(R,1); t_rawY = cell(R,1);
a_rawX = cell(R,1); a_rawY = cell(R,1);

for r=1:R
    A = Sall.(fn{r});
    if size(A,2)<2, error('%s: en az iki sÃ¼tun (t, ax) olmalÄ±.', fn{r}); end
    t0 = A(:,1); ax0 = A(:,2); ay0 = []; if size(A,2)>=3, ay0 = A(:,3); end
    [t0,iu] = unique(t0,'stable'); ax0 = ax0(iu); if ~isempty(ay0), ay0=ay0(iu); end
    
    [tX,ax] = regrid_to_target(t0, ax0, prep);
    if ~isempty(ay0), [~,ay] = regrid_to_target(t0, ay0, prep); else, ay=[]; end
    
    t_rawX{r} = tX; a_rawX{r} = ax;
    t_rawY{r} = tX; a_rawY{r} = ay;   % Y varsa aynÄ± Ä±zgara
    
    % Hedef dt kontrol uyarÄ±sÄ± (resample_mode='off' iken dahi)
    dtX = median(diff(tX),'omitnan');
    tol = max(prep.tol_rel*max(prep.target_dt,eps), 1e-12);
    if abs(dtX - prep.target_dt) > tol
        warning('KayÄ±t #%d dt=%.6g s, hedef=%.6g s (resample kapalÄ±).', r, dtX, prep.target_dt);
    end
end

%% -------------------- Åiddet eÅŸitleme (scaled set) --------------------

t_sclX = t_rawX; t_sclY = t_rawY;
a_sclX = a_rawX; a_sclY = a_rawY;   % default: scaled = raw

if pp.on.intensity
    zeta_SA = pp.PSA.zeta; band_fac = pp.PSA.band_fac; Np_band = pp.PSA.Np_band;
    
    % Opsiyonel CMS hedefi
    useCMS = false; T_cms = []; Sa_cms = [];
    if pp.on.CMS && exist('cms_target.mat','file')
        Scms = load('cms_target.mat');
        if isfield(Scms,'T_cms') && isfield(Scms,'Sa_cms') && numel(Scms.T_cms)==numel(Scms.Sa_cms)
            T_cms  = Scms.T_cms(:); Sa_cms = Scms.Sa_cms(:); useCMS = true;
        end
    end
    
    % Ã–nceden hesaplanan PSA var mÄ±?
    params_vec = [T1, zeta_SA, band_fac(:).', Np_band, ...
        pp.PSA.downsample_dt, double(pp.PSA.downsample_enabled)];
    prep_key = hash_psa_prep(t_rawX, a_rawX, params_vec);

    if ~isempty(env.PSA_PREP_CACHE) && isfield(env.PSA_PREP_CACHE,'key') ...
            && strcmp(env.PSA_PREP_CACHE.key, prep_key)
        tPSA   = env.PSA_PREP_CACHE.tPSA;
        agPSA  = env.PSA_PREP_CACHE.agPSA;
        Sa_band = env.PSA_PREP_CACHE.Sa_band;
    else
        % Band Sa (parfor opsiyonel)
        tPSA   = cell(R,1);
        agPSA  = cell(R,1);
        Sa_band = zeros(R,1);
        canPar  = pp.PSA.use_parfor && ~isempty(ver('parallel'));
        if canPar
            parfor r=1:R
                if pp.PSA.downsample_enabled
                    [tLoc,agLoc] = psa_grid(t_rawX{r}, a_rawX{r}, pp.PSA.downsample_dt);
                else
                    tLoc  = t_rawX{r};
                    agLoc = a_rawX{r};
                end
                Sa_band(r) = f_band(tLoc, agLoc, T1, zeta_SA, band_fac, Np_band);
                tPSA{r}  = tLoc;
                agPSA{r} = agLoc;
            end
        else
            for r=1:R
                if pp.PSA.downsample_enabled
                    [tLoc,agLoc] = psa_grid(t_rawX{r}, a_rawX{r}, pp.PSA.downsample_dt);
                else
                    tLoc  = t_rawX{r};
                    agLoc = a_rawX{r};
                end
                Sa_band(r) = f_band(tLoc, agLoc, T1, zeta_SA, band_fac, Np_band);
                tPSA{r}  = tLoc;
                agPSA{r} = agLoc;
            end
        end
        env.PSA_PREP_CACHE = struct('key',prep_key,'tPSA',{tPSA}, ...
            'agPSA',{agPSA},'Sa_band',Sa_band);
    end
    
    Sa_band_target = median(Sa_band);
    
    % Her kayÄ±t iÃ§in Ã¶lÃ§ek uygula (yalnÄ±z GA amaÃ§lÄ±; solver varsayÄ±lanÄ± RAW)
    for r=1:R
        Sab_r  = f_band(tPSA{r}, agPSA{r}, T1, zeta_SA, band_fac, Np_band);
        s_band = Sa_band_target / max(Sab_r,eps);
        
        if useCMS
            Sa_rec = f_vec(tPSA{r}, agPSA{r}, T_cms, zeta_SA);
            nume   = max(sum(Sa_rec.*Sa_cms),eps);
            deno   = max(sum(Sa_rec.*Sa_rec),eps);
            s_cms  = nume/deno;
            s_hyb  = s_band^(1-pp.gammaCMS) * s_cms^(pp.gammaCMS);
        else
            s_hyb  = s_band;
        end
        
        a_sclX{r} = s_hyb * a_rawX{r};
        if ~isempty(a_rawY{r}), a_sclY{r} = s_hyb * a_rawY{r}; end
    end
end

%% -------------------- Arias pencereleri (raw & scaled) -----------------

[t5x_raw,t95x_raw,t5y_raw,t95y_raw] = deal(zeros(R,1),zeros(R,1),nan(R,1),nan(R,1));
[t5x_scl,t95x_scl,t5y_scl,t95y_scl] = deal(zeros(R,1),zeros(R,1),nan(R,1),nan(R,1));

for r=1:R
    if pp.on.arias
        [t5x_raw(r), t95x_raw(r)] = arias_win(t_rawX{r}, a_rawX{r}, 0.05, 0.95);
        if ~isempty(a_rawY{r}), [t5y_raw(r), t95y_raw(r)] = arias_win(t_rawY{r}, a_rawY{r}, 0.05, 0.95); end
        
        [t5x_scl(r), t95x_scl(r)] = arias_win(t_sclX{r}, a_sclX{r}, 0.05, 0.95);
        if ~isempty(a_sclY{r}), [t5y_scl(r), t95y_scl(r)] = arias_win(t_sclY{r}, a_sclY{r}, 0.05, 0.95); end
    else
        t5x_raw(r)=t_rawX{r}(1);  t95x_raw(r)=t_rawX{r}(end);
        t5x_scl(r)=t_sclX{r}(1);  t95x_scl(r)=t_sclX{r}(end);
        if ~isempty(a_rawY{r}), t5y_raw(r)=t5x_raw(r);  t95y_raw(r)=t95x_raw(r); end
        if ~isempty(a_sclY{r}), t5y_scl(r)=t5x_scl(r);  t95y_scl(r)=t95x_scl(r); end
    end
end

%% -------------------- Baz parametre setleri (parametrebulur uyumlu) ---

% Damper geometri + malzeme (BAZ)
geom.Dp    = 0.125;                 % piston Ã§apÄ± [m]
geom.Lgap  = 0.0055;                 % etkin strok boÅŸluÄŸu [m]
geom.d_o   = 2.7e-3;             % orifis Ã§apÄ± [m] (kÃ¼Ã§Ã¼ltÃ¼lmÃ¼ÅŸ)
geom.Lori  = 0.1;                 % orifis uzunluÄŸu [m]
geom.Kd    = 1.6e8;                % yaÄŸÄ±n k_b (sÄ±kÄ±ÅŸabilirlikten gelen) iÃ§in Ã¶lÃ§ek [Pa]
geom.Ebody = 2.1e11;               % gÃ¶vde elastisite modÃ¼lÃ¼ [Pa]
geom.Ap    = pi*geom.Dp^2/4;       % piston alanÄ± [m^2] (tÃ¼retilen)

% Yay (spiral) malzeme/geo (BAZ)
sh.G      = 79e9;                  % kayma modÃ¼lÃ¼ [Pa]
sh.d_w    = 12e-4;                 % tel Ã§apÄ± [m]
sh.D_m    = 80e-2;                  % yay ortalama Ã§apÄ± [m]
sh.n_turn = 56;                    % sarÄ±m sayÄ±sÄ± [-]


% Orifis / akışkan / termal parametre structları
orf   = struct();
hyd   = struct();
therm = struct();
num   = struct();
if ~exist('cfg','var') || ~isstruct(cfg), cfg = struct(); end

% Numerik guardlar (parametrebulur ile aynı isimler)
num.dP_cap       = getfield_default(num,'dP_cap',       1e7);
num.mu_min_phys  = getfield_default(num,'mu_min_phys',  0.25);
num.softmin_eps  = getfield_default(num,'softmin_eps',  1e3);

% Damper parametrelerini türet
hyd.nStories = n - 1;
[geom, orf, hyd, therm] = init_damper_params(geom, orf, hyd, therm);
nd = hyd.n_parallel;

% Qsat/doygunluk için referans debi (ilk kodla uyumlu)
cd_ref        = max(orf.CdInf, orf.Cd0);
dp_for_qcap   = getfield_default(num,'dP_cap', 3e8);     % num.dP_cap varsa onu kullan
Ae_ref        = max(cd_ref * orf.Ao_eff, 1e-12);
num.Qcap_big  = getfield_default(num,'Qcap_big', ...
    hyd.Vmin_fac * Ae_ref * sqrt( 2*dp_for_qcap / max(therm.rho_ref,100) ) );

% Yay rijitlikleri
k_p   = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);   % coil spring
k_h   = geom.Kd*geom.Ap^2/geom.Lgap;            % hidrolik (tek damper)
k_s   = geom.Ebody*geom.Ap/geom.Lgap;           % gövde (tek damper)
k_hyd = 1/(1/max(k_h,eps) + 1/max(k_s,eps));    % seri birleşim
k_sd  = nd * (k_hyd + k_p);                     % kat başına efektif yay

% Bilgi: referans laminer direnç (T_ref'te)
R_lam0 = (128*therm.mu_ref*geom.Lori/(pi*geom.d_o^4)) / (nd * orf.n_orf);


%% -------------------- SeÃ§im yardÄ±mcÄ±larÄ± ------------------------------

pickA = @(SRC,Rid,DIR) pick_series(SRC,Rid,DIR, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl);

%% ==================== GA KoÅŸusu (isteÄŸe baÄŸlÄ±) ========================

% Ã‡ok aÅŸamalÄ± GA: Ã¶nce set-1 (geometri/rijitlik), sonra set-3 (akÄ±ÅŸkan/termal)
% KullanÄ±m: ga.opt.multi_stage  = [1 2 3 4 5]; // yoksa mevcut tek aÅŸamalÄ± akÄ±ÅŸ sÃ¼rer
if isfield(ga,'opt') && isfield(ga.opt,'enable') && ga.opt.enable
    % --- GA baÅŸlamadan: simulate ve yardÄ±mcÄ±lar iÅŸÃ§ilere eklensin ---
    assert(~isempty(which('simulate')), 'simulate.m path''te deÄŸil veya gÃ¶rÃ¼nmÃ¼yor.');
    pool = gcp('nocreate');
    if ~isempty(pool)
        addAttachedFiles(pool, { ...
            'simulate.m','ga_trial_log.m','trial_x_all28_from_structs.m', ...
            'decode_design_apply.m','ga_get_bounds.m','ensure_cfg_defaults.m', ...
            'mck_with_damper_adv.m','ga_call_compat.m' ...
            });
    end
    % --- AÅŸama listesi (yoksa tek aÅŸama: mevcut design_set) ---
    if isfield(ga.opt,'multi_stage') && ~isempty(ga.opt.multi_stage)
        stage_list = ga.opt.multi_stage(:).';
    else
        stage_list = ga.design_set;
    end
    
    last_best = []; last_set = NaN;
    
    for si = 1:numel(stage_list)
        ga.design_set = stage_list(si);
        
        
        % --- DeÄŸiÅŸken sÄ±nÄ±rlarÄ± + override ---
        [lb,ub,int_idx,~] = ga_get_bounds(ga.design_set);
        if isfield(ga.opt,'lb_override') && ~isempty(ga.opt.lb_override)
            lb = max(lb, ga.opt.lb_override);
        end
        
        if isfield(ga.opt,'ub_override') && ~isempty(ga.opt.ub_override)
            ub = min(ub, ga.opt.ub_override);
        end
        if isempty(int_idx), IntCon = [];
        else, IntCon = int_idx(:)'; end
        
        % --- Bilgi satÄ±rÄ± (banner) ---
        switch ga.design_set
            case 1
                fprintf(['[GA] set=1 | pop=%d | gen=%d | Cd0=%.2f..%.2f | CdInf=%.2f..%.2f | Rec=%g..%g | p_{exp}=%.2f..%.2f | cav_{sf}=%.2f..%.2f | K_{leak}=%.e..%.e | Vmin_{fac}=%.2f..%.2f | d_o=%.2f..%.2f mm | Lori=%.2f..%.2f mm\\n'], ...
                    ga.opt.popsize, ga.opt.maxgen, lb(1),ub(1), lb(2),ub(2), lb(3),ub(3), lb(4),ub(4), lb(5),ub(5), lb(6),ub(6), lb(7),ub(7), 1e3*lb(8),1e3*ub(8), 1e3*lb(9),1e3*ub(9));
            case 2
                fprintf(['[GA] set=2 | pop=%d | gen=%d | Dp=%.1f..%.1f mm | Lgap=%.1f..%.1f mm | Kd=%.2e..%.2e Pa | dP_{cap}=%.1e..%.1e Pa | F_{cap}=%.1e..%.1e N | n_{turn}=%d..%d | d_w=%.1f..%.1f mm | D_m=%.1f..%.1f mm | mu_{ref}=%.2f..%.2f Pa·s | beta_0=%.1e..%.1e Pa | resF=%g..%g\\n'], ...
                    ga.opt.popsize, ga.opt.maxgen, 1e3*lb(1),1e3*ub(1), 1e3*lb(2),1e3*ub(2), lb(3),ub(3), lb(4),ub(4), lb(5),ub(5), round(lb(6)),round(ub(6)), 1e3*lb(7),1e3*ub(7), 1e3*lb(8),1e3*ub(8), lb(9),ub(9), lb(10),ub(10), lb(11),ub(11));
            case 3
                fprintf(['[GA] set=3 | pop=%d | gen=%d | b_{\\mu}=%.4f..%.4f 1/°C | b_{\\beta}=%.4f..%.4f 1/°C | \\rho=%.0f..%.0f kg/m^3 | \\alpha_{\\rho}=%.1e..%.1e 1/°C | T_{ref}=%.0f..%.0f °C | Antoine(A,B,C)=(%.2f..%.2f, %.0f..%.0f, %.0f..%.0f)\\n'], ...
                    ga.opt.popsize, ga.opt.maxgen, lb(1),ub(1), lb(2),ub(2), lb(3),ub(3), lb(4),ub(4), lb(5),ub(5), lb(6),ub(6), lb(7),ub(7), lb(8),ub(8));
            case 4
                fprintf(['[GA] set=4 | pop=%d | gen=%d | hA_{os}=%g..%g W/K | hA_{o,env}=%g..%g | hA_{s,env}=%g..%g | c_{p,oil}=%g..%g | c_{p,steel}=%g..%g | T_{env}=%g..%g °C | T0=%g..%g °C | Ts0=%g..%g °C\\n'], ...
                    ga.opt.popsize, ga.opt.maxgen, lb(1),ub(1), lb(2),ub(2), lb(3),ub(3), lb(4),ub(4), lb(5),ub(5), lb(6),ub(6), lb(7),ub(7), lb(8),ub(8));
        end
        
        % --- BaÅŸlangÄ±Ã§ nÃ¼fusu (LHS) ---
        P0 = lhs_population(lb,ub,ga.opt.popsize);
        if isfield(ga.opt,'lhs_refine') && ga.opt.lhs_refine
            P0(1,:) = 0.5*(lb+ub);
        end
        
        % --- AÅŸama-Ã¶zel kÄ±sÄ±t/PF ayarÄ± ---
        cons_stage = cons;
        cfg_stage  = cfg;
        switch ga.design_set
            case 1  % Set1: Orifis+Akış
                cons_stage.on.dp_quant    = true;
                cons_stage.on.qsat_margin = true;
                cons_stage.on.cav_frac    = true;
                cons_stage.on.thermal_dT  = false;
                if ~isfield(cfg_stage,'solver') || ~isstruct(cfg_stage.solver), cfg_stage.solver = struct(); end
                cfg_stage.solver.preset = 'balanced';
                
            case 2  % Set2: Geometri+Hidrolik
                cons_stage.on.dp_quant    = true;
                cons_stage.on.qsat_margin = true;
                cons_stage.on.cav_frac    = true;
                cons_stage.on.thermal_dT  = false;
                if ~isfield(cfg_stage,'solver') || ~isstruct(cfg_stage.solver), cfg_stage.solver = struct(); end
                cfg_stage.solver.preset = 'balanced';
                
            case 3  % Set3: Termal Öz
                cons_stage.on.dp_quant    = true;
                cons_stage.on.qsat_margin = true;
                cons_stage.on.cav_frac    = true;
                cons_stage.on.thermal_dT  = true;
                if ~isfield(cfg_stage,'solver') || ~isstruct(cfg_stage.solver), cfg_stage.solver = struct(); end
                cfg_stage.solver.preset = 'tight';
                
            case 4  % Set4: Termal Sınır/IC + Isı geçişi
                cons_stage.on.dp_quant    = true;
                cons_stage.on.qsat_margin = true;
                cons_stage.on.cav_frac    = true;
                cons_stage.on.thermal_dT  = true;
                if ~isfield(cfg_stage,'solver') || ~isstruct(cfg_stage.solver), cfg_stage.solver = struct(); end
                cfg_stage.solver.preset = 'tight';
                
            case 5  % Set5: PF + Guards
                cons_stage.on.dp_quant    = true;
                cons_stage.on.qsat_margin = true;
                cons_stage.on.cav_frac    = true;
                cons_stage.on.thermal_dT  = false;
                if ~isfield(cfg_stage,'solver') || ~isstruct(cfg_stage.solver), cfg_stage.solver = struct(); end
                cfg_stage.solver.preset = 'balanced';
        end
        
        % --- Fitness sarÄ±cÄ± ---
        inner_fitfun = @(xx) eval_fitness_for_x(xx, ga.design_set, ...
            obj, cons_stage, pp.tail_sec, ...
            t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
            t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
            M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_stage, env, ...
            ga.opt.cache_enable, ga.opt.fail_early_k);
        
        
        
        fhandle = @(xx) compact_log_wrapper(xx, inner_fitfun);
        
        % --- GA OutputFcn init: nesil bazlÄ± en iyiyi kaydetmek iÃ§in ---
        save_path = '';
        if isfield(ga.opt,'save_log') && ~isempty(ga.opt.save_log)
            [pth,base,ext] = fileparts(ga.opt.save_log); if isempty(pth), pth='.'; end
            save_path = fullfile(pth, sprintf('%s_stage%d_current%s', base, ga.design_set, ext));
        end
        try, ga_log_outfun('init', struct('stage', ga.design_set, 'save_path', save_path)); end
        % Clear any pending soft-stop sentinel before starting GA
        try, ga_log_outfun('clearstop'); end
        % Also remove default STOP_GA file if left over
        try
            sf_def = fullfile('out','STOP_GA');
            if exist(sf_def,'file'), delete(sf_def); end
        end
        % Ensure trial log header exists once (client-side) to avoid worker truncation
        try, ga_trial_log('init_header', struct('file', fullfile('out','trial_log.csv'))); end
        
        % --- GA seÃ§enekleri ---
        opts = optimoptions('ga', ...
            'UseParallel', ga.opt.use_parallel, ...
            'InitialPopulationMatrix', P0, ...
            'PopulationSize', ga.opt.popsize, ...
            'MaxGenerations', ga.opt.maxgen, ...
            'Display','iter', ...
            'PlotFcn', {@gaplotbestf, @gaplotstopping, @gaplotrange, @gaplotbestindiv}, ...
            'PlotInterval', 1, ...
            'OutputFcn', @ga_log_outfun);
        
        
        % --- GÃ¼venli Ã¶n-tanÄ±mlar (xbest henÃ¼z yokken kullanmak YASAK) ---
        xbest = []; fbest = []; pop = []; scores = []; exitflag = []; output = struct();
        
        % --- Ã‡alÄ±ÅŸtÄ±r ---
        try
            [xbest, fbest, output, pop, scores, exitflag] = ga_call_compat(fhandle, lb, ub, IntCon, opts);
        catch ME
            % GA erken durdurulduysa, en iyi bireyi toplayÄ±p devam et
            fprintf('WARN: GA interrupted: %s\n', ME.message);
            best = [];
            try, best = ga_log_outfun('get'); end
            if ~isempty(best) && isfield(best,'x') && ~isempty(best.x)
                xbest = best.x(:).'; fbest = best.f; exitflag = NaN;
            else
                % PARETO Ã¼zerinden aynÄ± stage iÃ§in en iyiyi bul
                try
                    idx = find(env.PARETO.set == ga.design_set);
                    if ~isempty(idx)
                        [~,ik] = min(env.PARETO.F(idx));
                        xbest = env.PARETO.x{idx(ik)}; fbest = env.PARETO.F(idx(ik)); exitflag = NaN;
                    end
                end
            end
            if isempty(xbest)
                xbest = 0.5*(lb+ub); fbest = inf; exitflag = NaN;
            end
            output.message = 'interrupted'; pop = []; scores = [];
        end
        
        
        % --- Fallback: GA bir ÅŸey dÃ¶ndÃ¼rmediyse orta noktayÄ± kullan ---
        if isempty(xbest)
            xbest = 0.5*(lb+ub);
            fbest = inf;
            pop   = []; scores = [];
        end
        
        % --- Log & kalÄ±cÄ±laÅŸtÄ±r ---
        runLog = struct('stage',ga.design_set,'xbest',xbest,'fbest',fbest,'output',output, ...
            'pop',pop,'scores',scores,'exitflag',exitflag, ...
            'ga_options',opts,'bounds',struct('lb',lb,'ub',ub),'seed',ga.opt.seed);
        if isfield(ga.opt,'save_log') && ~isempty(ga.opt.save_log)
            try
                [pth,base,ext] = fileparts(ga.opt.save_log); if isempty(pth), pth='.'; end
                save(fullfile(pth, sprintf('%s_stage%d%s',base,ga.design_set,ext)), 'runLog');
            catch ME
                warning('runLog kaydedilemedi: %s', ME.message);
            end
        end
        
        % --- Bu aÅŸamanÄ±n en iyisini uygula â†’ bir sonraki aÅŸama bunun Ã¼stÃ¼nde arar ---
        if numel(xbest) ~= numel(lb)
            xbest = 0.5*(lb+ub);
        end
        ga.enable = true; ga.x = xbest;
        [geom, sh, orf, hyd, therm, num, ga, env] = decode_design_apply(ga, geom, sh, orf, hyd, therm, num, env);
        ga.enable = false; ga.x = [];
        
        last_best = xbest; last_set = ga.design_set;
        fprintf('\n=== GA Stage %d Bitti ===\nBest f = %.6g | xbest = [%s]\n', ...
            ga.design_set, fbest, join(string(xbest.'),', '));
        
        % --- SÄ±nÄ±r daraltma: bir SONRAKÄ° aÅŸama aynÄ± set ise
        if si < numel(stage_list)
            next_set = stage_list(si+1);
            try
                [lb_sh, ub_sh] = shrink_bounds_from_pop(pop, scores, lb, ub, ga.opt.keep_top, ga.opt.buffer_fac);
                if ~isempty(IntCon)
                    lb_sh(IntCon) = ceil(lb_sh(IntCon));
                    ub_sh(IntCon) = floor(ub_sh(IntCon));
                    lb_sh = max(lb, lb_sh);
                    ub_sh = min(ub, ub_sh);
                end
                if next_set == ga.design_set
                    ga.opt.lb_override = lb_sh;
                    ga.opt.ub_override = ub_sh;
                else
                    ga.opt.lb_override = [];
                    ga.opt.ub_override = [];
                end
            catch
                ga.opt.lb_override = [];
                ga.opt.ub_override = [];
            end
        end
        
        % --- AÅŸama bazlÄ± nesil sayÄ±sÄ± (varsa) ---
        if isfield(ga.opt,'maxgen_stage') && numel(ga.opt.maxgen_stage)>=si+1
            ga.opt.maxgen = ga.opt.maxgen_stage(si+1);
        end
    end
    
    
    % Son aÅŸamanÄ±n xâ€™i overlay/raporlar iÃ§in kalsÄ±n
    ga.enable = true; ga.x = last_best; ga.best_x = last_best; ga.best_set = last_set;
    ga_dbg = struct('enable',true,'design_set',ga.best_set,'x',ga.best_x);
    [geom_dbg, ~, orf_dbg, hyd_dbg, therm_dbg, ~, ~, env] = decode_design_apply(ga_dbg, geom, sh, orf, hyd, therm, num, env);
    
    fprintf('DBG set-%d: n_orf=%d | d_o=%.3f mm | Ao=%.3e m^2 | Lgap=%.1f mm | Vmin_fac=%.2f | Lh=%.3f mm | resFactor=%.0f\n', ...
        ga.best_set, orf_dbg.n_orf, 1e3*geom_dbg.d_o, orf_dbg.n_orf*(pi*geom_dbg.d_o^2/4), ...
        1e3*geom_dbg.Lgap, hyd_dbg.Vmin_fac, 1e3*hyd_dbg.Lh, therm_dbg.resFactor);
    
end

% --- GA sonucunu uygula ve tÃ¼retilenleri gÃ¼ncelle ---
if isfield(ga, 'best_x') && ~isempty(ga.best_x)
    ga.enable = true;
    ga.x = ga.best_x;
    ga.design_set = ga.best_set;
    [geom, sh, orf, hyd, therm, num, ga, env] = decode_design_apply(ga, geom, sh, orf, hyd, therm, num, env);
    ga.enable = false;
    ga.x = [];
end


% Türetilen parametreleri güncelle
hyd.nStories = n - 1;
[geom, orf, hyd, therm] = init_damper_params(geom, orf, hyd, therm);
nd = hyd.n_parallel;
cd_ref        = max(orf.CdInf, orf.Cd0);
dp_for_qcap   = getfield_default(num,'dP_cap', 3e8);
Ae_ref        = max(cd_ref * orf.Ao_eff, 1e-12);
num.Qcap_big  = getfield_default(num,'Qcap_big', ...
    hyd.Vmin_fac * Ae_ref * sqrt( 2*dp_for_qcap / max(therm.rho_ref,100) ) );
k_p   = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);
k_h   = geom.Kd*geom.Ap^2/geom.Lgap;
k_s   = geom.Ebody*geom.Ap/geom.Lgap;
k_hyd = 1/(1/max(k_h,eps) + 1/max(k_s,eps));
k_sd  = nd * (k_hyd + k_p);
R_lam0 = (128*therm.mu_ref*geom.Lori/(pi*geom.d_o^4)) / (nd * orf.n_orf);



%% -------------------- (1) Tek koÅŸu: gÃ¶rselleÅŸtirme/diagnostic ---------

rec = min(max(1, vis.sel_rec), R);
dir = upper(string(vis.sel_dir));
[t_sim, ag_sim, t5_sim, t95_sim] = pickA(sel.sim_source, rec, dir);
% Kuyruk ekle (sim)
dt_sim  = median(diff(t_sim));
t_tail  = (t_sim(end)+dt_sim:dt_sim:t_sim(end)+pp.tail_sec).';
t_sim   = [t_sim; t_tail];
ag_sim  = [ag_sim; zeros(size(t_tail))];

% PF rampasÄ± (sim penceresine gÃ¶re) â€” guard
cfg = set_pf_ton_if_auto(cfg, t5_sim, 0.5);

fprintf('SIM source=%s | rec #%d | dir=%s | N=%d | dt=%.3g s | Arias [%.2f, %.2f] s\n', ...
    sel.sim_source, rec, dir, numel(t_sim), dt_sim, t5_sim, t95_sim);

[x0_sim,a0_sim] = lin_MCK_consistent(t_sim, ag_sim, M, Cstr, K);
[xD_sim, aD_sim, dlog_sim, vD_sim] = mck_with_damper_adv( ...
    t_sim, ag_sim, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg);
% Grafik seti (opsiyonel ikinci koÅŸu)
use_plot_run = false;
cfg_plot = cfg;   % tern(...) Ã§aÄŸrÄ±sÄ± iÃ§in hazÄ±r dursun
if vis.dual_run_if_needed && ~strcmpi(sel.plot_source, sel.sim_source)
    [t_plot, ag_plot, t5_plot, t95_plot] = pickA(sel.plot_source, rec, dir);
    dt_plot = median(diff(t_plot));
    t_tail  = (t_plot(end)+dt_plot:dt_plot:t_plot(end)+pp.tail_sec).';
    t_plot  = [t_plot; t_tail];
    ag_plot = [ag_plot; zeros(size(t_tail))];
    
    cfg_plot = set_pf_ton_if_auto(cfg_plot, t5_plot, 0.5);
    
    [x0,a0] = lin_MCK_consistent(t_plot, ag_plot, M, Cstr, K);
    [xD,aD,dlog,vD] = mck_with_damper_adv( ...
        t_plot, ag_plot, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg_plot);
    
    use_plot_run = true;
    fprintf('PLOT source=%s | rec #%d | dir=%s | N=%d | dt=%.3g s | Arias [%.2f, %.2f] s\n', ...
        sel.plot_source, rec, dir, numel(t_plot), dt_plot, t5_plot, t95_plot);
else
    % Plot, sim ile aynÄ±
    t_plot=t_sim; ag_plot=ag_sim; t5_plot=t5_sim; t95_plot=t95_sim;
    x0=x0_sim; a0=a0_sim; xD=xD_sim; aD=aD_sim; dlog=dlog_sim; vD=vD_sim;
    dt_plot = dt_sim;   % baÅŸlÄ±ktaki dt iÃ§in
    
end

% ---- Enerji/diagnostik (sadece yazdÄ±rma) ----
active_cfg  = tern(use_plot_run, cfg_plot, cfg);
dlog_active = tern(use_plot_run, dlog, dlog_sim);

E_orif  = dlog_active.E_cum(end);                 % orifis enerjisi (J)
P_struc = sum( (vD * Cstr) .* vD, 2 );
E_struc = trapz(t_plot, P_struc);

fprintf('E_orifice = %.3e J | E_struct = %.3e J | oran = %.3f\n', ...
    E_orif, E_struc, E_orif/max(E_struc,eps));

P_mech = sum( (dlog_active.F_story .* (vD(:,2:end)-vD(:,1:end-1))), 2 );
fprintf('âŸ¨P_mechâŸ© = %.3e W (negatif ise net sÃ¶nÃ¼m)\n', mean(P_mech,'omitnan'));

fprintf('CHECK â†’ leak=%.2e, Lh=%.3e, Ao=%.2e m^2, PF=%s\n', ...
    hyd.K_leak, hyd.Lh, orf.n_orf*(pi*geom.d_o^2/4), active_cfg.PF.mode);

fprintf('dP_orf q95 = %.3e Pa | Q95 = %.3e m^3/s\n', ...
    prctile(dlog_active.dP_orf_time_max,95), dlog_active.Q_abs_p95);

if isfield(dlog_active,'cav_margin_min')
    fprintf('cav_margin_min = %.1f kPa\n', 1e-3*dlog_active.cav_margin_min);
end

% === Tek FIGÃœR (ilk koddaki stil): Ã¼stte yer deÄŸiÅŸtirme, altta mutlak ivme ===
if vis.make_plots
    j_disp = min(max(1, obj.idx_disp_story), size(xD,2));   % yer degistirme indeksi
    j_acc  = min(max(1, obj.idx_acc_story),  size(aD,2));   % ivme indeksi
    
    % Mutlak ivme (dampersiz/damperli)
    a_abs0 = a0 + ag_plot * ones(1, size(a0,2));
    a_absD = aD + ag_plot * ones(1, size(aD,2));
    
    figure('Color','w','Visible','on','Name', ...
        sprintf('%d. kat: dampersiz vs damperli', j_disp));
    
    % --- Ãœst: yer deÄŸiÅŸtirme ---
    subplot(2,1,1); hold on; grid on;
    plot(t_plot, x0(:, j_disp), 'k--', 'DisplayName','dampersiz');
    plot(t_plot, xD(:, j_disp), 'b-',  'DisplayName','damperli');
    ylabel(sprintf('x_{%d} [m]', j_disp)); legend('show','Location','best');
    title(sprintf('N=%d, dt=%.4fs | d_o=%.1f mm, gain=%.2f, \\tau=%.2f s', ...
        numel(t_plot), dt_plot, 1e3*geom.d_o, active_cfg.PF.gain, active_cfg.PF.tau));
    
    % --- Alt: mutlak ivme ---
    subplot(2,1,2); hold on; grid on;
    plot(t_plot, a_abs0(:, j_acc), 'k--', 'DisplayName','dampersiz');
    plot(t_plot, a_absD(:, j_acc), 'b-',  'DisplayName','damperli');
    ylabel(sprintf('a_{%d} [m/s^2]', j_acc)); xlabel('t [s]'); legend('show','Location','best');
    % Konsol Ã¶zeti: 3. kat ivme ve 10. kat yer deÄŸiÅŸtirme maksimumlarÄ±
    acc3_series = a_absD(:, j_acc);
    x10_series  = xD(:, j_disp);
    [acc3_max_abs, i_acc] = max(abs(acc3_series));
    [x10_max_abs,  i_x  ] = max(abs(x10_series));
    t_acc = t_plot(min(i_acc, numel(t_plot)));
    t_x   = t_plot(min(i_x,   numel(t_plot)));
    fprintf('PLOT METRICS: a_abs(j=%d) max = %.6g at t=%.3f s | x(j=%d) max = %.6g at t=%.3f s\n', ...
        j_acc, acc3_max_abs, t_acc, j_disp, x10_max_abs, t_x);
end

%% -------------------- (2) AmaÃ§ fonksiyonu â€” kompakt (ilk koddaki gibi) ----

% Girdi olarak ÅŸunlar zaten mevcut olmalÄ±:
% t_plot, dt_plot, vD, dlog_active, M, K, k_sd, obj, cons, orf, hyd, therm, num, t5_plot, t95_plot

% Kapasiteler
dPcap_eff = (isfield(num,'dP_cap') && isfinite(num.dP_cap)) * num.dP_cap + ...
    (~(isfield(num,'dP_cap') && isfinite(num.dP_cap))) * 3e8;

cd_ref   = max(orf.CdInf, orf.Cd0);
Ae_ref   = cd_ref * max(orf.Ao_eff, 1e-12);                 % paralel ve Cd dahil
rho_ref  = max(therm.rho_ref, 100);
Qcap_ref = hyd.Vmin_fac * Ae_ref * sqrt(2*dPcap_eff / rho_ref);

% GÃ¶zlenenler (dlog_active Ã¼zerinden)
Qp95_max = dlog_active.Q_abs_p95;                            % m^3/s
dp95_max = prctile(dlog_active.dP_orf_time_max, 95);         % Pa

% Uygunluk bayraklarÄ±
Qmargin = (isfield(cons,'hyd') && isfield(cons.hyd,'Q_margin') && isfinite(cons.hyd.Q_margin)) ...
    * cons.hyd.Q_margin + (~(isfield(cons,'hyd') && isfield(cons.hyd,'Q_margin'))) * 0.90;
okQ  = (Qp95_max <= Qmargin * Qcap_ref);
okdp = (dp95_max <= dPcap_eff);

% zeta_eq tahmini (1. mod)
[PHI,LAM] = eig(K,M);
WW = sqrt(builtin('diag', LAM));               % 'diag' fonksiyonu gÃ¶lgelenmesin diye builtin

w1   = min(WW);
phi1 = PHI(:, WW==w1); phi1 = phi1(:,1);
phi1 = phi1 / max(abs(phi1));
m_eff = phi1.'*M*phi1;
k_eff = phi1.'*K*phi1;

% Arias penceresinde izlenen katÄ±n hÄ±z RMS'i
j_mon = min(max(1, obj.idx_disp_story), size(vD,2));
maskA = (t_plot >= t5_plot) & (t_plot <= t95_plot);
v_rms = sqrt(mean(vD(maskA, j_mon).^2, 'omitnan'));
Xeq   = v_rms / max(w1, 1e-6);
Est   = 0.5 * (k_eff + k_sd) * Xeq^2;

% Mekanik gÃ¼Ã§ten efektif sÃ¶nÃ¼m (yalnÄ±z sÃ¶nÃ¼mleyici kÄ±sÄ±m)
P_mech  = sum( dlog_active.F_story .* (vD(:,2:end) - vD(:,1:end-1)), 2 );
P_diss  = max(-P_mech, 0);                                  % negatif â†’ sÃ¶nÃ¼m kabulÃ¼
zeta_eq = sum(P_diss) * dt_plot / max(4*pi*Est, eps);

% <P_mech> iÅŸareti (pozitif â‡’ net sÃ¶nÃ¼m)
Pmech_avg = -mean(P_mech, 'omitnan');
okP = (Pmech_avg > 0);

% Skor
score = zeta_eq ...
    - 0.7*double(~okP) ...
    - 0.5*double(~okQ) ...
    - 0.5*double(~okdp);

fprintf('\n=== AMAÃ‡ (kompakt) ===\n');
fprintf('zeta_eq=%.3f | <P_mech>_diss=%.3e W | Q95=%.3e | dp95=%.3e | ok=[Q %d, dp %d, P %d] | score=%.3f\n', ...
    zeta_eq, Pmech_avg, Qp95_max, dp95_max, okQ, okdp, okP, score);

% GA minimizasyonu ile uyum iÃ§in amaÃ§ deÄŸeri
J = -score;

%% -------------------- (3) KÄ±sÄ±tlar: Penalty hesap (R kayÄ±t) ----------
[Penalty, cons_detail] = evaluate_constraints_over_records( ...
    cons, cons.src_for_constraints, obj, pp.tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, tern(ga.enable,ga.design_set,0), tern(ga.enable,ga.x,[]));

Fitness = J + Penalty;
fprintf('Î»pen = [tau=%.2f, stroke=%.2f, dT=%.2f, cav=%.2f]\n', ...
    cons.pen.lambda.spring_tau, cons.pen.lambda.stroke, ...
    cons.pen.lambda.thermal_dT, cons.pen.lambda.cav_frac);
fprintf('PF: mode=%s, t_on=%.2fs, tau=%.2f, gain=%.2f\n', ...
    cfg.PF.mode, cfg.PF.t_on, cfg.PF.tau, cfg.PF.gain);



fprintf('\n================ KÄ±sÄ±t & Ceza SonuÃ§larÄ± ================\n');
fprintf('Penalty = %.6g  |  Fitness = J + Penalty = %.6g\n', Penalty, Fitness);
% === KÄ±sÄ±t/cihaz Ã¶zeti â†’ CSV ===
try, mkdir('out'); end

ratio = cons_detail.ratios;
okflag = @(r) (r<=1+1e-12);

Names   = {'spring_tau','spring_slender','stroke','force_cap','dp_quant','thermal_dT','cav_frac','qsat_margin'};
Limits  = [cons.spring.tau_allow, cons.spring.lambda_max, cons.stroke.util_factor*geom.Lgap, ...
    cons.force.F_cap, num.dP_cap, cons.thermal.cap_C, cons.hyd.cav_frac_cap, cons.hyd.Q_margin * getfield_default(num,'Qcap_big', ...
    0.4 * ( max(max(orf.CdInf,orf.Cd0)*orf.Ao_eff, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) )];

Ratios  = [ratio.spring_tau, ratio.spring_slender, ratio.stroke, ratio.force_cap, ...
    ratio.dp_quant, ratio.thermal_dT, ratio.cav_frac, ratio.qsat_margin];

Flags = strings(numel(Ratios),1);
okmask = Ratios <= (1 + 1e-12);
Flags(okmask)  = "OK";
Flags(~okmask) = "VIOL";


Vals_Fmax   = max(cons_detail.Fmax_records,   [], 'omitnan');
Vals_stroke = max(cons_detail.stroke_records, [], 'omitnan');
Vals_dpq    = max(cons_detail.dpq_records,    [], 'omitnan');
Vals_dT     = max(cons_detail.dT_records,     [], 'omitnan');
Vals_cav    = max(cons_detail.cav_records,    [], 'omitnan');
Vals_Qp95   = max(cons_detail.Qp95_records,   [], 'omitnan');

T = table(Names.', Ratios.', Flags, ...
    'VariableNames', {'Constraint','Ratio','Flag'});

% Ek cihaz metrikleri ayrÄ± tablo: (zarf deÄŸerleri)
T_dev = table(Vals_Fmax, Vals_stroke, Vals_dpq, Vals_dT, Vals_cav, Vals_Qp95, ...
    'VariableNames', {'Fmax_N','stroke_m','dp_q95_Pa','dT_C','cav95','Q95_m3s'});

writetable(T,    fullfile('out','cons_summary.csv'));
writetable(T_dev,fullfile('out','device_summary.csv'));
fprintf('CSV yazÄ±ldÄ±: out/cons_summary.csv, out/device_summary.csv\n');


% ---------- YardÄ±mcÄ±lar ----------
hinge  = @(r) max(0, r - 1);
norm0  = @(x) (abs(x) < 1e-12) * 0 + (abs(x) >= 1e-12) .* x;  % -0.000 yerine 0.000
pwr    = cons.pen.power;
lam    = cons.pen.lambda;

pen_sum = 0;   % bileÅŸenlerden yeniden hesaplanan toplam (kontrol amaÃ§lÄ±)

% ---------- Ã–zet oranlar + bireysel ceza katkÄ±larÄ± ----------
cdt = cons_detail;  % kÄ±saltma

if cons.on.spring_tau
    r = norm0(cdt.ratios.spring_tau);
    pen_i = lam.spring_tau * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('Ï„_max/Ï„_allow         = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.spring_slender
    r = norm0(cdt.ratios.spring_slender);
    pen_i = lam.spring_slender * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('L_free/D_m / Î»_max    = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.stroke
    r = norm0(cdt.ratios.stroke);
    pen_i = lam.stroke * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('stroke/(0.9*L_gap)    = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.force_cap
    r = norm0(cdt.ratios.force_cap);
    pen_i = lam.force_cap * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('F_max/F_cap           = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.dp_quant
    r = norm0(cdt.ratios.dp_quant);
    pen_i = lam.dp_quant * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('q_Î”p/dP_cap           = %.3f (q=%.3f)   [pen=%.3g]\n', r, cons.dp.q, pen_i);
end

if cons.on.thermal_dT
    r = norm0(cdt.ratios.thermal_dT);
    pen_i = lam.thermal_dT * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('Î”T/Î”T_cap             = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.cav_frac
    r = norm0(cdt.ratios.cav_frac);
    pen_i = lam.cav_frac * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('cav95/Î³_cap           = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.qsat_margin
    r = norm0(cdt.ratios.qsat_margin);
    pen_i = lam.qsat_margin * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('Q95/(margin*Qcap_big) = %.3f   [pen=%.3g]\n', r, pen_i);
end


% BigM (sim baÅŸarÄ±sÄ±zlÄ±ÄŸÄ±) bilgilendirmesi ve katkÄ±sÄ±
if cons.on.fail_bigM && cdt.any_fail
    pen_BigM = cons.pen.bigM * lam.fail_bigM;
    pen_sum  = pen_sum + pen_BigM;
    fprintf('Sim baÅŸarÄ±sÄ±z (â‰¥1 koÅŸu): BigM cezasÄ± eklendi. [pen=%.3g]\n', pen_BigM);
end

fprintf('--- Penalty (bileÅŸenlerden) â‰ˆ %.6g\n', pen_sum);



% ---------- Ä°steÄŸe baÄŸlÄ± tanÄ±lama: zarf/metrik Ã¶zetleri ----------
if isfield(cons_detail,'Fmax_records')
    % KayÄ±t zarfÄ± / agregeler (dpq_all iÃ§in agresyon kuralÄ±nÄ± uygula)
    Fmax_all   = max(cons_detail.Fmax_records,   [], 'omitnan');
    stroke_all = max(cons_detail.stroke_records, [], 'omitnan');
    if isfield(cons,'dp') && isfield(cons.dp,'agg') && strcmpi(cons.dp.agg,'cvar')
        dpq_all = cvar_from_samples(cons_detail.dpq_records(:), cons.alpha_CVaR_cons);
    else
        dpq_all = max(cons_detail.dpq_records, [], 'omitnan');
    end
    dT_all   = max(cons_detail.dT_records,   [], 'omitnan');
    cav_all  = max(cons_detail.cav_records,  [], 'omitnan');
    Qp95_all = max(cons_detail.Qp95_records, [], 'omitnan');
    
    fprintf('--- TanÄ±lama (zarf deÄŸerleri) ---\n');
    fprintf('Fmax=%.3e N | stroke=%.3e m | dpq=%.3e Pa | Î”T=%.2f C | cav95=%.3f | Q95=%.3e m^3/s\n', ...
        Fmax_all, stroke_all, dpq_all, dT_all, cav_all, Qp95_all);
end


% ======================================================================
%                             FONKSÄ°YONLAR
% ======================================================================
function [lb2,ub2] = shrink_bounds_from_pop(pop, scores, lb, ub, keep_top, buf)
if isempty(pop) || isempty(scores)
    lb2 = lb; ub2 = ub; return;
end
[~,ix] = sort(scores(:),'ascend');                 % en iyi kÃ¼Ã§Ã¼k
K = max(1, ceil(keep_top * size(pop,1)));
P = pop(ix(1:K), :);                               % elitler
p10 = prctile(P,10,1);
p90 = prctile(P,90,1);
span = max(p90 - p10, 1e-12);
lb2 = max(lb, p10 - buf.*span);
ub2 = min(ub, p90 + buf.*span);
end

function y = tern(cond,a,b)
if cond
    y = a;
else
    y = b;
end
end


function [J, out] = compute_objective_over_records( ...
    src, obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, varargin)

if numel(varargin) >= 2
    design_set = varargin{1}; x_ga = varargin{2};
else
    design_set = 0; x_ga = [];
end

if ~isempty(x_ga)
    try
        ga_local = struct('enable', true, 'design_set', design_set, 'x', x_ga);
        [~,~,~,~,~,num_ga,~, env] = decode_design_apply(ga_local, geom, sh, orf, hyd, therm, num, env);
        if isfield(num_ga, 'F_cap')
            cons.force.F_cap = num_ga.F_cap;
        end
    catch
    end
end
...
    % local_design_ratios_one_dir iÃ§inde:
% SimÃ¼lasyon penceresi: aynÄ± veriye kuyruk ekleyelim (PF ramp guard uyumu)
% Kaynak seÃ§ici
useScaled = strcmpi(src,'scaled');
if useScaled
    tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
    t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
else
    tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
    t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
end

R = numel(tX);
out_rec(R) = struct('d_rel',NaN,'a_rel',NaN,'J_r',NaN);
logs = {};
fail_ratio = 5;  % fallback multiplier when simulation fails

% Ã–nce REFERANSLAR: (damper yok, Î¼=1.0) â†’ her kayÄ±t & (X,Y mevcutsa) yÃ¶n iÃ§in sabitlenir
ref = struct('X',struct('d',nan(R,1),'a',nan(R,1)), ...
    'Y',struct('d',nan(R,1),'a',nan(R,1)));

for r=1:R
    % X yÃ¶nÃ¼ referans
    if ~isempty(aX{r})
        [dref, aref] = local_ref_metrics(tX{r}, aX{r}, t5x(r), t95x(r), obj, M,Cstr,K, n);
        ref.X.d(r)=dref; ref.X.a(r)=aref;
    end
    % Y yÃ¶nÃ¼ referans (varsa)
    if ~isempty(aY{r})
        [dref, aref] = local_ref_metrics(tY{r}, aY{r}, t5y(r), t95y(r), obj, M,Cstr,K, n);
        ref.Y.d(r)=dref; ref.Y.a(r)=aref;
    end
end

% Sonra TASARIM: her kayÄ±t â†’ yÃ¶n zarfÄ± â†’ Î¼-agg â†’ J_r
for r=1:R
    [d_rel_X, a_rel_X] = local_design_ratios_one_dir('X', r);
    [d_rel_Y, a_rel_Y] = local_design_ratios_one_dir('Y', r);
    
    switch lower(obj.dir_mode)
        case 'xonly', d_rel = d_rel_X; a_rel = a_rel_X;
        case 'yonly', d_rel = d_rel_Y; a_rel = a_rel_Y;
        otherwise     % 'envelope'
            d_rel = max([d_rel_X, d_rel_Y], [], 'omitnan');
            a_rel = max([a_rel_X, a_rel_Y], [], 'omitnan');
    end
    
    % AÄŸÄ±rlÄ±klÄ± toplam
    J_r = obj.weights_da(1)*d_rel + obj.weights_da(2)*a_rel;
    
    out_rec(r).d_rel = d_rel;
    out_rec(r).a_rel = a_rel;
    out_rec(r).J_r   = J_r;
end

% CVaR(Î±) hesap
alpha = min(max(obj.alpha_CVaR,eps),0.99);
Jlist = [out_rec.J_r].';
J = cvar_from_samples(Jlist, alpha);
if ~isfinite(J)
    msg = 'Objective J is NaN/Inf';
    logs{end+1} = msg; log_msg('error', msg);
    J = 1e9;
end

out = struct('records', out_rec, 'logs', {logs});

% ---- iÃ§ yardÄ±mcÄ±lar ----

    function [d_ref, a_ref] = local_ref_metrics(t, ag, t5, t95, obj, M,C,K, n)
        % Kuyruk eklemeden baz referans (damper yok)
        [x0,a0] = lin_MCK_consistent(t, ag, M, C, K);
        d_ref = max(abs(x0(:,min(obj.idx_disp_story,n))));
        a_ref = acc_metric_from_series(t, a0(:,min(obj.idx_acc_story,n)), t5, t95, obj);
        d_ref = max(d_ref, eps);
        a_ref = max(a_ref, eps);
    end

    function [d_rel_dir, a_rel_dir] = local_design_ratios_one_dir(which, r)
        % Skip unused direction when dir_mode is X-only/Y-only
        try
            dm = lower(obj.dir_mode);
        catch
            dm = 'xonly';
        end
        if (strcmpi(which,'Y') && strcmp(dm,'xonly')) || (strcmpi(which,'X') && strcmp(dm,'yonly'))
            d_rel_dir = NaN; a_rel_dir = NaN; return;
        end
        
        switch upper(string(which))
            case "X"
                t=tX{r}; ag=aX{r}; t5=t5x(r); t95=t95x(r);
                d_ref=ref.X.d(r); a_ref=ref.X.a(r);
            case "Y"
                t=tY{r}; ag=aY{r}; t5=t5y(r); t95=t95y(r);
                d_ref=ref.Y.d(r); a_ref=ref.Y.a(r);
        end
        if isempty(ag) || ~isfinite(d_ref) || ~isfinite(a_ref)
            d_rel_dir = NaN; a_rel_dir = NaN; return;
        end
        
        
        
        % --- SimÃ¼lasyon penceresi: aynÄ± veriye kuyruk ekleyelim ---
        
        dt = median(diff(t));
        tail_sec_loc = tail_sec;
        t_tail = (t(end)+dt : dt : t(end)+tail_sec_loc).';
        t_s    = [t;  t_tail];
        ag_s   = [ag; zeros(size(t_tail))];
        
        % PF ramp t_on: Arias t5 + 3
        cfg_dir = cfg;  % simulate will handle PF.t_on if auto_t_on is true
        
        
        % Î¼ senaryolarÄ±
        mus   = obj.mu_scenarios(:).';
        d_vals = nan(size(mus));  a_vals = nan(size(mus));
        
        for k = 1:numel(mus)
            [resp, env] = simulate( ...
                design_set, x_ga, mus(k), t_s, ag_s, ...
                M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_dir, env);
            
            if ~resp.ok
                msg = sprintf('simulate fail (rec=%d dir=%s mu=%g): %s', r, which, mus(k), resp.msg);
                logs{end+1} = msg; %#ok<AGROW>
                log_msg('warn', msg);
                d_vals(k) = fail_ratio * d_ref;
                a_vals(k) = fail_ratio * a_ref;
                continue;
            end
            
            x  = resp.y.x(1:numel(t), :);   % kuyruÄŸu at
            aR = resp.y.a(1:numel(t), :);
            
            d  = max(abs(x(:,min(obj.idx_disp_story,n))));
            am = acc_metric_from_series(t, aR(:,min(obj.idx_acc_story,n)), t5, t95, obj);
            
            d_vals(k) = d;  a_vals(k) = am;
        end
        
        % Î¼-aggregation
        switch lower(obj.mu_aggregate)
            case 'weighted'
                w = obj.mu_weights(:); w = w/sum(w);
                d_agg = nansum(w.*d_vals(:));
                a_agg = nansum(w.*a_vals(:));
            otherwise
                d_agg = max(d_vals, [], 'omitnan');
                a_agg = max(a_vals, [], 'omitnan');
        end
        
        d_rel_dir = d_agg / max(d_ref, eps);
        a_rel_dir = a_agg / max(a_ref, eps);
        if ~isfinite(d_rel_dir)
            d_rel_dir = fail_ratio;
            msg = sprintf('NaN d_rel at rec=%d dir=%s', r, which);
            logs{end+1} = msg; log_msg('warn', msg);
        end
        if ~isfinite(a_rel_dir)
            a_rel_dir = fail_ratio;
            msg = sprintf('NaN a_rel at rec=%d dir=%s', r, which);
            logs{end+1} = msg; log_msg('warn', msg);
        end
    end

end

function val = acc_metric_from_series(t, a, t5, t95, obj)
% Arias penceresi iÃ§i metrikler
if obj.use_arias_window
    w = (t>=t5 & t<=t95);
else
    w = true(size(t));
end
aa = a(w); tt = t(w);
if isempty(aa) || numel(aa)<2
    val = NaN; return;
end
switch lower(obj.acc_metric)
    case 'energy'
        val = trapz(tt, aa.^2);                 % enerji
    case 'rms+p95'
        rmsv = sqrt(mean(aa.^2,'omitnan'));
        p95  = prctile(abs(aa),95);
        val  = rmsv + obj.p95_penalty_w * p95;  % hibrit
    otherwise % 'rms'
        val = sqrt(mean(aa.^2,'omitnan'));
end
end

function v = cvar_from_samples(x, alpha)
% Ã–rneklerden CVaR_Î± (Average Value-at-Risk)
x = x(:); x = x(isfinite(x));
if isempty(x), v = NaN; return; end
q = quantile(x, 1-alpha);
tail = x(x>=q);  % bÃ¼yÃ¼k-kÃ¶tÃ¼ kuyruk
if isempty(tail), v = q; else, v = mean(tail); end
end
function [J1, out] = compute_J1_IDR_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, varargin)

% Optional GA args + precomputed references (ref_pre)
ref_pre = [];
if numel(varargin) >= 3
    design_set = varargin{1}; x_ga = varargin{2}; ref_pre = varargin{3};
elseif numel(varargin) >= 2
    design_set = varargin{1}; x_ga = varargin{2};
else
    design_set = 0; x_ga = [];
end

% Kaynak seÃ§imi
useScaled = strcmpi(src,'scaled');
if useScaled
    tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
    t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
else
    tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
    t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
end

R  = numel(tX);
Ns = n - 1;
if isscalar(h_story_m), h_story = repmat(h_story_m, Ns, 1);
else,                   h_story = h_story_m(:); end

out = struct('d_rel', nan(R,1));
mus = obj.mu_scenarios(:).';
isWeighted = strcmpi(obj.mu_aggregate,'weighted');
if isWeighted
    wmu = obj.mu_weights(:); wmu = wmu / max(sum(wmu), eps);
end

for r = 1:R
    % --- X yÃ¶nÃ¼ referans IDR
    dref_X = NaN; dagg_X = NaN;
    if ~isempty(aX{r}) && ~strcmpi(obj.dir_mode,'yonly')
        [dref_X, dagg_X] = local_idr_ref_and_agg(tX{r}, aX{r}, t5x(r), t95x(r), 'X', r);
    end
    
    % --- Y yÃ¶nÃ¼ (varsa)
    dref_Y = NaN; dagg_Y = NaN;
    if ~isempty(aY{r}) && ~strcmpi(obj.dir_mode,'xonly')
        [dref_Y, dagg_Y] = local_idr_ref_and_agg(tY{r}, aY{r}, t5y(r), t95y(r), 'Y', r);
    end
    
    % --- YÃ¶n zarfÄ±
    switch lower(obj.dir_mode)
        case 'xonly', d_ref = dref_X; d_agg = dagg_X;
        case 'yonly', d_ref = dref_Y; d_agg = dagg_Y;
        otherwise
            d_ref = max([dref_X, dref_Y], [], 'omitnan');
            d_agg = max([dagg_X, dagg_Y], [], 'omitnan');
    end
    
    out.d_rel(r) = d_agg / max(d_ref, eps);
end

J1 = cvar_from_samples(out.d_rel, obj.alpha_CVaR);

% ---- local helper
    function [d_ref, d_agg] = local_idr_ref_and_agg(t, ag, t5, t95, dirStr, rid)
        % Referans (damper yok, Î¼=1.0), pencere iÃ§inde IDR zarfÄ±
        w = (t >= t5 & t <= t95);
        if ~isempty(ref_pre)
            if strcmpi(dirStr,'X'), d_ref = ref_pre.X.d(rid); else, d_ref = ref_pre.Y.d(rid); end
        else
            [x0,~] = lin_MCK_consistent(t, ag, M, Cstr, K);
            drift0 = x0(:,2:end) - x0(:,1:end-1);                       % Nt x Ns
            idr0   = abs(drift0(w,:)) ./ (ones(sum(w),1) * h_story(:)'); % Nt_w x Ns
            d_ref  = max(idr0,[],'all');                                 % skaler
        end
        
        % TasarÄ±m: kuyruk ekle + PF ramp guard
        dt     = median(diff(t));
        t_tail = (t(end)+dt : dt : t(end)+tail_sec).';
        t_s    = [t; t_tail];
        ag_s   = [ag; zeros(size(t_tail))];
        
        cfg_dir = cfg;  % simulate will handle PF.t_on if auto_t_on is true
        
        vals = nan(size(mus));
        for k = 1:numel(mus)
            [resp, env] = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir, env);
            if ~resp.ok
                vals(k) = 5 * d_ref;  % fail durumunda gÃ¼venli bÃ¼yÃ¼k ceza
                continue;
            end
            xD   = resp.y.x(1:numel(t), :);     % kuyruÄŸu at
            driftD = xD(:,2:end) - xD(:,1:end-1);
            idrD   = abs(driftD(w,:)) ./ (ones(sum(w),1) * h_story(:)');
            vals(k)= max(idrD,[],'all');
        end
        
        if isWeighted, d_agg = nansum(wmu .* vals(:));
        else,          d_agg = max(vals, [], 'omitnan'); end
    end
end


function [J2, out] = compute_J2_ACC_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, varargin)

ref_pre = [];
if numel(varargin) >= 3
    design_set = varargin{1}; x_ga = varargin{2}; ref_pre = varargin{3};
elseif numel(varargin) >= 2
    design_set = varargin{1}; x_ga = varargin{2};
else, design_set = 0; x_ga = []; end


useScaled = strcmpi(src,'scaled');
if useScaled
    tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
    t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
else
    tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
    t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
end

R   = numel(tX);
w_p = obj.p95_penalty_w;
out = struct('a_rel', nan(R,1));

mus = obj.mu_scenarios(:).';
isWeighted = strcmpi(obj.mu_aggregate,'weighted');
if isWeighted
    wmu = obj.mu_weights(:); wmu = wmu / max(sum(wmu), eps);
end

for r = 1:R
    % --- X yÃ¶nÃ¼: A_env ref ve agg
    aref_X = NaN; aagg_X = NaN;
    if ~isempty(aX{r}) && ~strcmpi(obj.dir_mode,'yonly')
        [aref_X, aagg_X] = local_acc_ref_and_agg(tX{r}, aX{r}, t5x(r), t95x(r), 'X', r);
    end
    
    % --- Y yÃ¶nÃ¼ (varsa)
    aref_Y = NaN; aagg_Y = NaN;
    if ~isempty(aY{r}) && ~strcmpi(obj.dir_mode,'xonly')
        [aref_Y, aagg_Y] = local_acc_ref_and_agg(tY{r}, aY{r}, t5y(r), t95y(r), 'Y', r);
    end
    
    % --- YÃ¶n zarfÄ±
    switch lower(obj.dir_mode)
        case 'xonly', a_ref = aref_X; a_agg = aagg_X;
        case 'yonly', a_ref = aref_Y; a_agg = aagg_Y;
        otherwise
            a_ref = max([aref_X, aref_Y], [], 'omitnan');
            a_agg = max([aagg_X, aagg_Y], [], 'omitnan');
    end
    
    out.a_rel(r) = a_agg / max(a_ref, eps);
end

J2 = cvar_from_samples(out.a_rel, obj.alpha_CVaR);

% ---- local helper
    function [A_ref, A_agg] = local_acc_ref_and_agg(t, ag, t5, t95, dirStr, rid)
        % Referans (damper yok) mutlak ivme zarfÄ± (RMS+p95)
        w = (t >= t5 & t <= t95);
        if ~isempty(ref_pre)
            if strcmpi(dirStr,'X'), A_ref = ref_pre.X.a(rid); else, A_ref = ref_pre.Y.a(rid); end
        else
            [~,a_rel0] = lin_MCK_consistent(t, ag, M, Cstr, K);  % relatif
            a_abs0 = a_rel0 + ag(:) * ones(1, size(a_rel0,2));  % her kat iÃ§in mutlak
            A_i = zeros(1,size(a_abs0,2));
            for i=1:size(a_abs0,2)
                ai = a_abs0(w,i);
                rmsv = sqrt(mean(ai.^2,'omitnan'));
                p95  = prctile(abs(ai),95);
                A_i(i) = rmsv + w_p * p95;
            end
            A_ref = max(A_i);   % kat zarfÄ±
        end
        
        % TasarÄ±m: kuyruk ekle + PF guard
        dt     = median(diff(t));
        t_tail = (t(end)+dt : dt : t(end)+tail_sec).';
        t_s    = [t; t_tail];
        ag_s   = [ag; zeros(size(t_tail))];
        cfg_dir = cfg;  % simulate will handle PF.t_on if auto_t_on is true
        
        vals = nan(size(mus));
        for k = 1:numel(mus)
            [resp, env] = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir, env);
            if ~resp.ok
                vals(k) = 5 * A_ref;  % fail durumunda gÃ¼venli bÃ¼yÃ¼k ceza
                continue;
            end
            a_relD = resp.y.a(1:numel(t), :);           % relatif, kuyruÄŸu at
            a_absD = a_relD + ag(:) * ones(1,size(a_relD,2));
            
            Ai = zeros(1,size(a_absD,2));
            for i=1:size(a_absD,2)
                ai = a_absD(w,i);
                rmsv = sqrt(mean(ai.^2,'omitnan'));
                p95  = prctile(abs(ai),95);
                Ai(i)= rmsv + w_p * p95;
            end
            vals(k) = max(Ai);   % kat zarfÄ±
        end
        
        if isWeighted, A_agg = nansum(wmu .* vals(:));
        else,          A_agg = max(vals, [], 'omitnan'); end
    end
end

function [Penalty, out] = evaluate_constraints_over_records( ...
    cons, src, obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, varargin)

if numel(varargin) >= 2
    design_set = varargin{1}; x_ga = varargin{2};
else
    design_set = 0; x_ga = [];
end
...
    
% Kaynak seÃ§imi
useScaled = strcmpi(src,'scaled');
if useScaled
    tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
    t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
else
    tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
    t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
end

R   = numel(tX);
mus = cons.mu_scenarios(:).';
logs = {};

% ----- Erken Ã§Ä±kÄ±ÅŸ sayaÃ§ ayarlarÄ± -----
fail_early_k = 0;
if isfield(cons,'fail_early_k') && ~isempty(cons.fail_early_k) && isfinite(cons.fail_early_k)
    fail_early_k = max(0, round(cons.fail_early_k));
end
fail_count_global = 0;
early_break_all   = false;

% ToplayÄ±cÄ±lar
Fmax_records   = nan(R,1);
stroke_records = nan(R,1);
dpq_records    = nan(R,1);
dT_records     = nan(R,1);
cav_records    = nan(R,1);
Qp95_records   = nan(R,1);
any_fail_rec   = false(R,1);

for r = 1:R
    % ---- YÃ¶n zarfÄ± iÃ§in X/Y Ã¶lÃ§Ã¼leri
    rec_metrics = struct('Fmax',[],'stroke',[],'dpq',[],'dT',[],'cav95',[],'Qp95',[],'fail',[]);
    
    % Active directions: honor obj.dir_mode to avoid unused branch
    if strcmpi(obj.dir_mode,'xonly')
        dirs = ["X"];
    elseif strcmpi(obj.dir_mode,'yonly')
        dirs = ["Y"];
    else
        dirs = ["X","Y"];
    end
    for dir = dirs
        % Y yÃ¶nÃ¼ yoksa atla
        if dir=="Y" && (isempty(aY{r}) || all(isnan(aY{r})))
            continue;
        end
        
        % KayÄ±t/yÃ¶n serileri
        if dir=="X"
            t = tX{r}; ag = aX{r}; t5=t5x(r); t95=t95x(r);
        else
            t = tY{r}; ag = aY{r}; t5=t5y(r); t95=t95y(r);
        end
        
        % ---- Tail ekle (simÃ¼lasyon penceresi)
        dt    = median(diff(t));
        t_tail= (t(end)+dt : dt : t(end)+tail_sec).';
        t_s   = [t; t_tail];
        ag_s  = [ag; zeros(size(t_tail))];
        
        % ---- PF ramp korumasÄ±
        cfg_dir = cfg;  % simulate will handle PF.t_on if auto_t_on is true
        
        
        % ---- Î¼-senaryosu zarfÄ±
        Fmax_mu   = -Inf; stroke_mu = -Inf; dpq_mu = -Inf; dT_mu = -Inf; cav_mu = -Inf; Qp95_mu = -Inf;
        fail_mu   = false;
        
        % >>>>>>>>>>>>>>>>> Î¼ DÃ–NGÃœSÃœ (BURADA) <<<<<<<<<<<<<<<<<
        for k = 1:numel(mus)
            [resp, env] = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir, env);
            
            if ~resp.ok
                fail_mu = true;
                fail_count_global = fail_count_global + 1;
                msg = sprintf('simulate fail (rec=%d dir=%s mu=%g): %s', r, dir, mus(k), resp.msg);
                logs{end+1} = msg; %#ok<AGROW>
                log_msg('warn', msg);
                
                if fail_early_k > 0 && fail_count_global >= fail_early_k
                    early_break_all = true;   % tÃ¼m dÃ¶ngÃ¼lerden Ã§Ä±k
                    break;                    % Î¼-dÃ¶ngÃ¼sÃ¼nÃ¼ kÄ±r
                end
                continue;  % sÄ±radaki Î¼
            end
            
            % ---- Ã¶lÃ§Ã¼mler (baÅŸarÄ±lÄ± koÅŸu)
            Fmax_mu   = max(Fmax_mu,   resp.F_max);
            stroke_mu = max(stroke_mu, resp.stroke_max);
            dpq_mu    = max(dpq_mu,    resp.dP_q_time(cons.dp.q));
            dT_mu     = max(dT_mu,     resp.dT_est);
            mask = (t_s >= t5) & (t_s <= t95);
            if any(mask)
                cav_p95_win = prctile(resp.diag.cav_frac_t(mask), 95);
            else
                cav_p95_win = prctile(resp.diag.cav_frac_t, 95);  % emniyetli geri dÃ¶nÃ¼ÅŸ
            end
            cav_mu = max(cav_mu, cav_p95_win);
            Qp95_mu   = max(Qp95_mu,   resp.metrics.Q_abs_p95);
        end
        
        % >>>>>>>>>>>>>>>>> Î¼ DÃ–NGÃœSÃœ BÄ°TTÄ° <<<<<<<<<<<<<<<<<<<<
        
        % Î¼-dÃ¶ngÃ¼sÃ¼ eÅŸik nedeniyle kÄ±rÄ±ldÄ±ysa, yÃ¶n dÃ¶ngÃ¼sÃ¼nÃ¼ de bÄ±rak
        if early_break_all
            break;
        end
        
        % Bu yÃ¶nÃ¼n metriklerini biriktir
        rec_metrics.Fmax   = [rec_metrics.Fmax,   Fmax_mu];
        rec_metrics.stroke = [rec_metrics.stroke, stroke_mu];
        rec_metrics.dpq    = [rec_metrics.dpq,    dpq_mu];
        rec_metrics.dT     = [rec_metrics.dT,     dT_mu];
        rec_metrics.cav95  = [rec_metrics.cav95,  cav_mu];
        rec_metrics.Qp95   = [rec_metrics.Qp95,   Qp95_mu];
        rec_metrics.fail   = [rec_metrics.fail,   fail_mu];
    end
    
    % YÃ¶n dÃ¶ngÃ¼sÃ¼nden erken Ã§Ä±ktÄ±ysak kalan kayÄ±tlarÄ± INF kabul et
    if early_break_all
        any_fail_rec(r:end) = true;
        break;
    end
    
    if early_break_all
        Penalty = cons.pen.bigM * cons.pen.lambda.fail_bigM;
        out = struct();
        out.ratios = struct('spring_tau',0,'spring_slender',0,'stroke',0, ...
            'force_cap',0,'dp_quant',0,'thermal_dT',0, ...
            'cav_frac',0,'qsat_margin',0);
        out.any_fail = true;
        out.dpq_records    = dpq_records;
        out.Fmax_records   = Fmax_records;
        out.stroke_records = stroke_records;
        out.dT_records     = dT_records;
        out.cav_records    = cav_records;
        out.Qp95_records   = Qp95_records;
        out.logs           = logs;
        return
    end
    
    
    % ---- YÃ¶n zarfÄ± (X,Y â†’ max)
    Fmax_records(r)   = max(rec_metrics.Fmax,   [], 'omitnan');
    stroke_records(r) = max(rec_metrics.stroke, [], 'omitnan');
    dpq_records(r)    = max(rec_metrics.dpq,    [], 'omitnan');
    dT_records(r)     = max(rec_metrics.dT,     [], 'omitnan');
    cav_records(r)    = max(rec_metrics.cav95,  [], 'omitnan');
    Qp95_records(r)   = max(rec_metrics.Qp95,   [], 'omitnan');
    any_fail_rec(r)   = any(rec_metrics.fail);
end


% ---- KayÄ±t agregasyonlarÄ±
% Fmax ve stroke en-kÃ¶tÃ¼ kayÄ±t
Fmax_all   = max(Fmax_records,   [], 'omitnan');
stroke_all = max(stroke_records, [], 'omitnan');
dT_all     = max(dT_records,     [], 'omitnan');
cav_all    = max(cav_records,    [], 'omitnan');
Qp95_all   = max(Qp95_records,   [], 'omitnan');

% Î”p quantile: kayÄ±tlar arasÄ± 'max' veya 'CVaR'
switch lower(cons.dp.agg)
    case 'cvar'
        dpq_all = cvar_from_samples(dpq_records(:), cons.alpha_CVaR_cons);
    otherwise
        dpq_all = max(dpq_records, [], 'omitnan');
end

agg_vals = [Fmax_all stroke_all dpq_all dT_all cav_all Qp95_all];
if any(~isfinite(agg_vals))
    names = {'Fmax','stroke','dpq','dT','cav','Qp95'};
    msg = ['NaN/Inf aggregate: ' strjoin(names(~isfinite(agg_vals)), ', ')];
    logs{end+1} = msg; %#ok<AGROW>
    log_msg('error', msg);
    Penalty = cons.pen.bigM;
    out = struct('ratios',struct(),'any_fail',true, ...
        'dpq_records',dpq_records,'Fmax_records',Fmax_records, ...
        'stroke_records',stroke_records,'dT_records',dT_records, ...
        'cav_records',cav_records,'Qp95_records',Qp95_records, ...
        'logs',{logs});
    return;
end

% ---- Oranlar (â‰¥1 â†’ ihlal)
ratios = struct();

% K1: yay kesme gerilmesi (yalnÄ±z yay kolu kuvveti)
if cons.on.spring_tau
    k_p_est = sh.G*sh.d_w^4 / (8*sh.n_turn*sh.D_m^3);           % [N/m]
    Fspring_max = k_p_est * stroke_all;                          % [N]
    ratios.spring_tau = spring_tau_ratio(Fspring_max, sh, cons.spring.tau_allow);
else
    ratios.spring_tau = 0;
end

% K2: burkulma/serbest boy oranÄ±
if cons.on.spring_slender
    if cons.spring.use_fixed_length && isfinite(cons.spring.L_free_fixed)
        L_free = cons.spring.L_free_fixed;
    else
        L_free = cons.spring.L_free_auto_fac * geom.Lgap; % yaklaÅŸÄ±klama
    end
    lambda = (L_free / max(sh.D_m,eps)) / max(cons.spring.lambda_max,eps);
    ratios.spring_slender = lambda;
else
    ratios.spring_slender = 0;
end

% K3: strok
if cons.on.stroke
    ratios.stroke = stroke_all / max(cons.stroke.util_factor*geom.Lgap, eps);
else
    ratios.stroke = 0;
end

% K4: cihaz kuvvet limiti
if cons.on.force_cap
    ratios.force_cap = Fmax_all / max(cons.force.F_cap, eps);
else
    ratios.force_cap = 0;
end

% K5: Î”p quantile
if cons.on.dp_quant
    ratios.dp_quant = dpq_all / max(num.dP_cap, eps);
else
    ratios.dp_quant = 0;
end

% K6: termal Î”T
if cons.on.thermal_dT
    ratios.thermal_dT = dT_all / max(cons.thermal.cap_C, eps);
else
    ratios.thermal_dT = 0;
end

% K7: kavitasyon
if cons.on.cav_frac
    ratios.cav_frac = cav_all / max(cons.hyd.cav_frac_cap, eps);
else
    ratios.cav_frac = 0;
end

% K8: Q satÃ¼rasyon marjÄ±
if cons.on.qsat_margin
    Qcap = getfield_default(num,'Qcap_big', ...
        0.4 * ( max(max(orf.CdInf,orf.Cd0)*orf.Ao_eff, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
    ratios.qsat_margin = Qp95_all / max(cons.hyd.Q_margin * Qcap, eps);
    
else
    ratios.qsat_margin = 0;
end

% ---- Ceza hesabÄ±
hinge = @(r) max(0, r - 1);
pwr   = cons.pen.power;

pen = 0;
if cons.on.spring_tau,     pen = pen + cons.pen.lambda.spring_tau     * hinge(ratios.spring_tau    )^pwr; end
if cons.on.spring_slender, pen = pen + cons.pen.lambda.spring_slender * hinge(ratios.spring_slender)^pwr; end
if cons.on.stroke,         pen = pen + cons.pen.lambda.stroke         * hinge(ratios.stroke        )^pwr; end
if cons.on.force_cap
    pen = pen + cons.pen.lambda.force_cap      * hinge(ratios.force_cap     )^pwr; end
if cons.on.dp_quant,       pen = pen + cons.pen.lambda.dp_quant       * hinge(ratios.dp_quant      )^pwr; end
if cons.on.thermal_dT,     pen = pen + cons.pen.lambda.thermal_dT     * hinge(ratios.thermal_dT    )^pwr; end
if cons.on.cav_frac,       pen = pen + cons.pen.lambda.cav_frac       * hinge(ratios.cav_frac      )^pwr; end
if cons.on.qsat_margin,    pen = pen + cons.pen.lambda.qsat_margin    * hinge(ratios.qsat_margin   )^pwr; end

any_fail = any(any_fail_rec);
if cons.on.fail_bigM && any_fail
    pen = pen + cons.pen.bigM * cons.pen.lambda.fail_bigM;
end

Penalty = pen;

% Ã‡Ä±kÄ±ÅŸ detaylarÄ±
out = struct();
out.ratios    = ratios;
out.any_fail  = any_fail;
out.dpq_records = dpq_records;
out.Fmax_records = Fmax_records;
out.stroke_records = stroke_records;
out.dT_records = dT_records;
out.cav_records = cav_records;
out.Qp95_records = Qp95_records;
out.logs = logs;
end
function ratio = spring_tau_ratio(Fmax, sh, tau_allow)
C  = max(sh.D_m, eps) / max(sh.d_w, eps);
Kw = (4*C - 1)/(4*C - 4) + 0.615/C;
tau_max = (8 * Fmax * max(sh.D_m,eps) / max(pi*sh.d_w^3, eps)) * Kw;
ratio = tau_max / max(eps, tau_allow);
end

function P = lhs_population(lb,ub,N)
D = numel(lb); P = zeros(N,D);
for d=1:D
    edges = linspace(0,1,N+1);
    centers = (edges(1:end-1)+edges(2:end))/2;
    P(:,d) = lb(d) + centers(randperm(N))' .* (ub(d)-lb(d));
end
end

function [f, aux] = eval_fitness_for_x(x, design_set, ...
    obj, cons, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, ...
    use_cache, fail_early_k)

% ---- Ã–nbellek anahtarÄ±
persistent memo
if isempty(memo), memo = containers.Map('KeyType','char','ValueType','any'); end
key = sprintf('set%d|%s', design_set, mat2str(x,8));

if use_cache && isKey(memo,key)
    data = memo(key); f = data.f; aux = data.aux; return;
end

% ---- TasarÄ±ma uygula ve deÄŸerlendir
try
    % GA decodeâ€™u simulate iÃ§inde deÄŸil, doÄŸrudan burada yapmaya gerek yok;
    % compute/evaluate fonksiyonlarÄ± simulateâ€™i Ã§aÄŸÄ±rÄ±rken set/xâ€™Ä± geÃ§iriyoruz.
    
    % AmaÃ§: J
    goal_src = tern(obj.use_scaled_for_goal,'scaled','raw');
    [J, ~] = compute_objective_over_records( ...
        goal_src, obj, tail_sec, ...
        t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
        t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
        M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
        env, design_set, x);
    
    % KÄ±sÄ±t: Penalty (erken Ã§Ä±kÄ±ÅŸ desteÄŸi ile)
    cons_loc = cons; cons_loc.fail_early_k = fail_early_k;
    [Penalty, cons_detail] = evaluate_constraints_over_records( ...
        cons_loc, cons.src_for_constraints, obj, tail_sec, ...
        t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
        t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
        M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
        env, design_set, x);
    
    % === J1 & J2 (split) hesap â€” Pareto gÃ¼nlÃ¼ÄŸÃ¼ iÃ§in ===
    [J1_split, J2_split] = compute_objectives_split( ...
        tern(obj.use_scaled_for_goal,'scaled','raw'), obj, tail_sec, ...
        t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
        t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
        M,Cstr,K,n,geom,sh,orf,hyd,therm,num, ...
        cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
        env, design_set, x);
    
    % aux iÃ§ine koy
    aux_J1 = J1_split; aux_J2 = J2_split;
    
    % --- Pareto gÃ¼nlÃ¼ÄŸÃ¼ne yaz ---
    if ~isfield(env,'PARETO') || ~isstruct(env.PARETO)
        env.PARETO = struct('J1',[],'J2',[],'F',[],'Pen',[],...
            'set',[],'x',{{}},'feas',[]);
    end
    env.PARETO.J1(end+1,1)  = aux_J1;
    env.PARETO.J2(end+1,1)  = aux_J2;
    env.PARETO.F(end+1,1)   = J + Penalty;
    env.PARETO.Pen(end+1,1) = Penalty;
    env.PARETO.set(end+1,1) = design_set;
    env.PARETO.x{end+1,1}   = x(:).';
    env.PARETO.feas(end+1,1)= (Penalty <= 1e-6);    % eÅŸik: cezasÄ±z â‰ˆ fizibÄ±l
    
    f = J + Penalty;
    aux = struct('J',J,'Penalty',Penalty,'cons',cons_detail);
    
    % ---- Trial log (row per evaluation) ----
    try
        % 1) Map full 28-length design snapshot from currently used structs
        try
            ga_local=struct('enable',true,'design_set',design_set,'x',x);
            [geom2,sh2,orf2,hyd2,therm2,num2,~, env]=decode_design_apply(ga_local, geom, sh, orf, hyd, therm, num, env);
        catch
            geom2=geom; sh2=sh; orf2=orf; hyd2=hyd; therm2=therm; num2=num;
        end
        % PF overrides (if any) applied via cfg_with_ga for logging
        cfg_log = cfg_with_ga(design_set, x, geom2, sh2, orf2, hyd2, therm2, num2, cfg);
        vars32 = struct( ...
            'cav_sf',     getfield_default(orf2,'cav_sf',NaN), ...
            'resFactor',  getfield_default(therm2,'resFactor',NaN), ...
            'Lgap',       getfield_default(geom2,'Lgap',NaN), ...
            'n_turn',     getfield_default(sh2,'n_turn',NaN), ...
            'Cd0',        getfield_default(orf2,'Cd0',NaN), ...
            'p_exp',      getfield_default(orf2,'p_exp',NaN), ...
            'Lh',         getfield_default(hyd2,'Lh',NaN), ...
            'CdInf',      getfield_default(orf2,'CdInf',NaN), ...
            'K_leak',     getfield_default(hyd2,'K_leak',NaN), ...
            'mu_ref',     getfield_default(therm2,'mu_ref',NaN), ...
            'Vmin_fac',   getfield_default(hyd2,'Vmin_fac',NaN), ...
            'hA_os',      getfield_default(therm2,'hA_os',NaN), ...
            'Rec',        getfield_default(orf2,'Rec',NaN), ...
            'Kd',         getfield_default(geom2,'Kd',NaN), ...
            'beta0',      getfield_default(therm2,'beta0',NaN), ...
            'Dp',         getfield_default(geom2,'Dp',NaN), ...
            'b_mu',       getfield_default(therm2,'b_mu',NaN), ...
            'dP_cap',     getfield_default(num2,'dP_cap',NaN), ...
            'F_cap',      getfield_default(num2,'F_cap',NaN), ...
            'antoine_A',  getfield_default(therm2,'antoine_A',NaN), ...
            'antoine_B',  getfield_default(therm2,'antoine_B',NaN), ...
            'antoine_C',  getfield_default(therm2,'antoine_C',NaN), ...
            'T0_C',       getfield_default(therm2,'T0_C',NaN), ...
            'Ts0_C',      getfield_default(therm2,'Ts0_C',NaN), ...
            'T_env_C',    getfield_default(therm2,'T_env_C',NaN), ...
            'hA_o_env',   getfield_default(therm2,'hA_o_env',NaN), ...
            'hA_s_env',   getfield_default(therm2,'hA_s_env',NaN), ...
            'cp_oil',     getfield_default(therm2,'cp_oil',NaN), ...
            'cp_steel',   getfield_default(therm2,'cp_steel',NaN), ...
            'rho_ref',    getfield_default(therm2,'rho_ref',NaN), ...
            'T_ref_C',    getfield_default(therm2,'T_ref_C',NaN), ...
            'alpha_rho',  getfield_default(therm2,'alpha_rho',NaN), ...
            'b_beta',     getfield_default(therm2,'b_beta',NaN), ...
            'pf_tau',     getfield_default(cfg_log.PF,'tau',NaN), ...
            'pf_gain',    getfield_default(cfg_log.PF,'gain',NaN), ...
            'pf_t_on',    getfield_default(cfg_log.PF,'t_on',NaN), ...
            'p_amb',      getfield_default(orf2,'p_amb',NaN), ...
            'mu_min_phys',getfield_default(num2,'mu_min_phys',NaN), ...
            'softmin_eps',getfield_default(num2,'softmin_eps',NaN) ...
            );
        % 2) Penalty components
        ratios = cons_detail.ratios; pwr = cons.pen.power; lam = cons.pen.lambda;
        pen = struct();
        if cons.on.spring_tau,     pen.spring_tau     = lam.spring_tau     * max(0, ratios.spring_tau    -1)^pwr; else, pen.spring_tau=0; end
        if cons.on.spring_slender, pen.spring_slender = lam.spring_slender * max(0, ratios.spring_slender-1)^pwr; else, pen.spring_slender=0; end
        if cons.on.stroke,         pen.stroke         = lam.stroke         * max(0, ratios.stroke        -1)^pwr; else, pen.stroke=0; end
        if cons.on.force_cap,      pen.force_cap      = lam.force_cap      * max(0, ratios.force_cap     -1)^pwr; else, pen.force_cap=0; end
        if cons.on.dp_quant,       pen.dp_quant       = lam.dp_quant       * max(0, ratios.dp_quant      -1)^pwr; else, pen.dp_quant=0; end
        if cons.on.thermal_dT,     pen.thermal_dT     = lam.thermal_dT     * max(0, ratios.thermal_dT    -1)^pwr; else, pen.thermal_dT=0; end
        if cons.on.cav_frac,       pen.cav_frac       = lam.cav_frac       * max(0, ratios.cav_frac      -1)^pwr; else, pen.cav_frac=0; end
        if cons.on.qsat_margin,    pen.qsat_margin    = lam.qsat_margin    * max(0, ratios.qsat_margin   -1)^pwr; else, pen.qsat_margin=0; end
        pen.fail_bigM = (cons.on.fail_bigM && cons_detail.any_fail) * (cons.pen.bigM * lam.fail_bigM);
        
        % 3) zeta_eq, okP (single representative run: first X record, mu=1.0)
        % choose data source consistent with goal
        useScaled = strcmpi(tern(obj.use_scaled_for_goal,'scaled','raw'),'scaled');
        if useScaled
            t_cells = t_sclX; a_cells = a_sclX; t5s = t5x_scl; t95s = t95x_scl;
        else
            t_cells = t_rawX; a_cells = a_rawX; t5s = t5x_raw; t95s = t95x_raw;
        end
        zeta_eq = NaN; okP = NaN; okQ = ratios.qsat_margin <= 1+1e-12; okDp = ratios.dp_quant <= 1+1e-12;
        acc3_max_abs = NaN; x10_max_abs = NaN; acc3_area_pos=NaN; acc3_area_neg_abs=NaN; x10_area_pos=NaN; x10_area_neg_abs=NaN;
        rpick = 1;
        if ~isempty(t_cells) && numel(t_cells)>=rpick && ~isempty(a_cells{rpick})
            t0 = t_cells{rpick}; ag0 = a_cells{rpick}; t5 = t5s(rpick); t95 = t95s(rpick);
            % extend tail briefly for PF consistency
            dt = median(diff(t0)); t_tail = (t0(end)+dt:dt:t0(end)+tail_sec).';
            t_s = [t0; t_tail]; ag_s = [ag0; zeros(size(t_tail))];
            [resp, env] = simulate(design_set, x, 1.0, t_s, ag_s, M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg, env);
            if resp.ok
                % mode data
                [PHI,LAM] = eig(K,M); WW = sqrt(builtin('diag',LAM)); w1 = min(WW);
                phi1 = PHI(:, WW==w1); phi1 = phi1(:,1); phi1 = phi1 / max(abs(phi1));
                m_eff = phi1.'*M*phi1; k_eff = phi1.'*K*phi1;
                % k_sd estimate
                nd = max(1, getfield_default(hyd,'n_parallel',1));
                Ap = pi*geom.Dp^2/4;
                k_p = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);
                k_h = geom.Kd*Ap^2/geom.Lgap; k_s = geom.Ebody*Ap/geom.Lgap;
                k_hyd = 1/(1/max(k_h,eps) + 1/max(k_s,eps));
                k_sd  = nd * (k_hyd + k_p);
                % build signals without tail
                Nt = numel(t0);
                vD = resp.y.v(1:Nt,:);
                dvel = vD(:,2:end) - vD(:,1:end-1);
                F_story = resp.diag.F_story(1:Nt,:);
                P_mech = sum(F_story .* dvel, 2);
                P_diss = max(-P_mech, 0);
                v_rms = sqrt(mean(vD(t0>=t5 & t0<=t95, min(max(1,obj.idx_disp_story), size(vD,2))).^2, 'omitnan'));
                Xeq = v_rms / max(w1, 1e-6);
                Est = 0.5 * (k_eff + k_sd) * Xeq^2;
                zeta_eq = sum(P_diss) * median(diff(t0)) / max(4*pi*Est, eps);
                okP = (-mean(P_mech,'omitnan')) > 0;
                % absolute acceleration and displacement metrics
                a_rel = resp.y.a(1:Nt,:); a_abs = a_rel + ag0(:)*ones(1,size(a_rel,2));
                j3  = min(3, size(a_abs,2)); j10 = min(10, size(resp.y.x,2));
                acc3 = a_abs(:, j3);
                x10  = resp.y.x(1:Nt, j10);
                acc3_max_abs = max(abs(acc3)); x10_max_abs = max(abs(x10));
                dt0 = median(diff(t0));
                acc3_pos = max(acc3, 0); acc3_neg = min(acc3, 0);
                x10_pos  = max(x10, 0);  x10_neg  = min(x10, 0);
                acc3_area_pos     = trapz(t0, acc3_pos);
                acc3_area_neg_abs = trapz(t0, -acc3_neg);
                x10_area_pos      = trapz(t0, x10_pos);
                x10_area_neg_abs  = trapz(t0, -x10_neg);
            end
        end
        
        score_val = zeta_eq - 0.7*double(~okP) - 0.5*double(~okQ) - 0.5*double(~okDp);
        ga_trial_log('append', struct( ...
            'design_set', design_set, 'ga_x', x(:).', 'vars32', vars32, ...
            'zeta_eq', zeta_eq, 'okP', okP, 'okQ', okQ, 'okDp', okDp, ...
            'Penalty', Penalty, 'J', J, 'J1', aux_J1, 'J2', aux_J2, 'score', score_val, 'pen', pen, ...
            'acc3_max_abs', acc3_max_abs, 'x10_max_abs', x10_max_abs, ...
            'acc3_area_pos', acc3_area_pos, 'acc3_area_neg_abs', acc3_area_neg_abs, ...
            'x10_area_pos', x10_area_pos, 'x10_area_neg_abs', x10_area_neg_abs));
    catch
        % logging should not break optimization
    end
    
    if use_cache, memo(key) = struct('f',f,'aux',aux); end
catch ME
    % GÃ¼venli bÃ¼yÃ¼k ceza + ayrÄ±ntÄ±lÄ± rapor
    f = 1e9;
    aux = struct('err', getReport(ME, 'extended', 'hyperlinks', 'off'));
end

end

function cfg2 = cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg)
%#ok<INUSD>
cfg2 = ensure_cfg_defaults(cfg);
try
    if ~isempty(x) && isfinite(design_set) && design_set==5
        % x = [PF_tau, PF_gain, PF_t_on, p_amb, mu_min_phys, softmin_eps]
        if ~isfield(cfg2,'PF') || ~isstruct(cfg2.PF), cfg2.PF = struct(); end
        cfg2.PF.tau = x(1);
        cfg2.PF.gain = x(2);
        cfg2.PF.t_on = x(3);
    end
catch
end
end

% =============================== Alt YapÄ± ===============================
function f = compact_log_wrapper(x, inner_fitfun)
% Tek satÄ±r log: [idx  #feval  J  (J+Penalty)  nViol]
% inner_fitfun: [f, aux] dÃ¶ndÃ¼rebilir (f=J+Penalty), aux.J ve aux.cons.ratios iÃ§erebilir.
persistent ROW FEVAL
if isempty(ROW),   ROW   = 0; end
if isempty(FEVAL), FEVAL = 0; end
ROW   = ROW + 1;
FEVAL = FEVAL + 1;

J = NaN; nviol = 0;

try
    [f_val, aux] = inner_fitfun(x);   % iki Ã§Ä±ktÄ± destekli
    if isstruct(aux)
        if isfield(aux,'J'), J = aux.J; end
        if isfield(aux,'cons') && isfield(aux.cons,'ratios') && isstruct(aux.cons.ratios)
            fn = fieldnames(aux.cons.ratios);
            vals = zeros(numel(fn),1);
            for i=1:numel(fn), vals(i) = aux.cons.ratios.(fn{i}); end
            nviol = sum(vals > 1+1e-12);
        end
    end
    f = f_val;
    
    % Hata raporu geldiyse konsola bas
    if isstruct(aux) && isfield(aux,'err') && ~isempty(aux.err)
        fprintf('ERR: %s\n', aux.err);
    end
    
catch
    % inner_fitfun tek Ã§Ä±ktÄ± verirse
    try
        f = inner_fitfun(x);
    catch
        f = 1e9;   % gÃ¼venli bÃ¼yÃ¼k ceza
    end
end

fprintf('%6d %15d %13.3f %13.3f %8d\n', ROW, FEVAL, J, f, nviol);
end


function [t2,a2] = regrid_to_target(t1,a1,prep)
% Tekil zaman dÃ¼zelt, hedef dt'ye gÃ¶re (auto/off/force)
[t1,iu]=unique(t1,'stable'); a1=a1(iu);
dt1 = median(diff(t1),'omitnan');
tgt = prep.target_dt; tol = max(prep.tol_rel*max(tgt,eps), 1e-12);
switch lower(prep.resample_mode)
    case 'off'
        t2=t1; a2=a1;
    case 'force'
        t2 = (t1(1):tgt:t1(end)).';
        a2 = interp1(t1,a1,t2,prep.regrid_method,'extrap');
    otherwise % 'auto'
        if abs(dt1 - tgt) <= tol
            t2=t1; a2=a1;
        else
            t2 = (t1(1):tgt:t1(end)).';
            a2 = interp1(t1,a1,t2,prep.regrid_method,'extrap');
            warning('Resample: dt=%.6gâ†’%.6g s | N=%d', dt1, tgt, numel(t2));
        end
end
end

function [t,a,t5,t95] = pick_series(src, rid, dir, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl)

src = lower(string(src)); dir = upper(string(dir));
switch src
    case "raw"
        if dir=="X", t=t_rawX{rid}; a=a_rawX{rid}; t5=t5x_raw(rid); t95=t95x_raw(rid);
        else,         t=t_rawY{rid}; a=a_rawY{rid}; t5=t5y_raw(rid); t95=t95y_raw(rid); end
    otherwise % 'scaled'
        if dir=="X", t=t_sclX{rid}; a=a_sclX{rid}; t5=t5x_scl(rid); t95=t95x_scl(rid);
        else,         t=t_sclY{rid}; a=a_sclY{rid}; t5=t5y_scl(rid); t95=t95y_scl(rid); end
end
if isempty(a), error('KayÄ±t #%d iÃ§in %s yÃ¶nÃ¼ mevcut deÄŸil.', rid, dir); end
end

function [M,K,C] = make_KCM(n,mv,kv,cv)
M = diag(mv); K = zeros(n); C = zeros(n);
for i=1:n
    kL=kv(i); cL=cv(i);
    if i<n, kU=kv(i+1); cU=cv(i+1); else, kU=0; cU=0; end
    K(i,i) = kL + (i<n)*kU;   C(i,i) = cL + (i<n)*cU;
    if i>1, K(i,i-1)=-kL; C(i,i-1)=-cL; end
    if i<n, K(i,i+1)=-kU; C(i,i+1)=-cU; end
end
end

function [x,a_rel] = lin_MCK_consistent(t, ag, M, C, K)
% Persistent cache for baseline linear MCK response
persistent memo
if isempty(memo)
    memo = containers.Map('KeyType','char','ValueType','any');
end

n  = size(M,1);
% Build a robust key from sizes/sums and the input series
sig = [size(M,1); size(M,2); sum(M(:)); sum(C(:)); sum(K(:))];
vec = [double(n); sig(:); t(:); ag(:)];
try
    key = md5(num2str(vec(:)', '%.16g,'));
catch
    key = sprintf('n%d|N%d|t%.6g..%.6g|ag%.6g..%.6g|sM%.6g|sC%.6g|sK%.6g', ...
        n, numel(t), t(1), t(end), ag(1), ag(end), sum(M(:)), sum(C(:)), sum(K(:)));
end

if isKey(memo, key)
    data = memo(key);
    x = data.x; a_rel = data.a_rel; return;
end

r  = ones(n,1);
dt = median(diff(t));
agf = griddedInterpolant(t,ag,'linear','nearest');
odef = @(tt,z) [ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - M*r*agf(tt) ) ];
z0 = zeros(2*n,1);
opts = odeset('RelTol',2e-3,'AbsTol',1e-6,'MaxStep',max(dt*10,2e-3),'InitialStep',max(dt*0.25,1e-3));
sol = ode23tb(odef,[t(1) t(end)],z0,opts);
t_end = sol.x(end); idx = find(t <= t_end + 1e-12);
if isempty(idx), x=nan(numel(t),n); a_rel=x; warning('lin_MCK_consistent: early stop'); return; end
t_use = t(idx); Z = deval(sol,t_use).';
x_use = Z(:,1:n); v_use = Z(:,n+1:end);
a_use = ( -(M\(C*v_use.' + K*x_use.')).' - ag(1:numel(t_use)).*r.' );
x=nan(numel(t),n); a_rel=x; x(1:numel(t_use),:)=x_use; a_rel(1:numel(t_use),:)=a_use;

memo(key) = struct('x',x,'a_rel',a_rel);
end



function [J1, J2] = compute_objectives_split( ...
    src, obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, design_set, x_ga)

% Compute objective components separately:
%   J1 - IDR zarf oranÄ± CVaR
%   J2 - mutlak ivme (RMS+p95) zarf oranÄ± CVaR
h_story_m = 3.0 * ones(n-1,1); % zaten Ã¼stte de var; buraya da koyduk

% Precompute references once per record/direction and pass down
useScaled = strcmpi(src,'scaled');
if useScaled
    tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
    t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
else
    tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
    t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
end

R  = numel(tX);
if isscalar(h_story_m), h_story = repmat(h_story_m, n-1, 1); else, h_story = h_story_m(:); end
w_p = obj.p95_penalty_w;

refIDR = struct('X',struct('d',nan(R,1)),'Y',struct('d',nan(R,1)));
refACC = struct('X',struct('a',nan(R,1)),'Y',struct('a',nan(R,1)));

for r = 1:R
    % X direction
    if ~isempty(aX{r})
        [x0, a_rel0] = lin_MCK_consistent(tX{r}, aX{r}, M, Cstr, K);
        w = (tX{r} >= t5x(r) & tX{r} <= t95x(r));
        % IDR ref
        drift0 = x0(:,2:end) - x0(:,1:end-1);
        idr0   = abs(drift0(w,:)) ./ (ones(sum(w),1) * h_story(:)');
        refIDR.X.d(r) = max(idr0,[],'all');
        % ACC ref
        a_abs0 = a_rel0 + aX{r}(:) * ones(1, size(a_rel0,2));
        A_i = zeros(1, size(a_abs0,2));
        for i=1:size(a_abs0,2)
            ai = a_abs0(w,i);
            rmsv = sqrt(mean(ai.^2,'omitnan'));
            p95  = prctile(abs(ai),95);
            A_i(i) = rmsv + w_p * p95;
        end
        refACC.X.a(r) = max(A_i);
    end
    % Y direction (if present)
    if ~isempty(aY{r})
        [x0, a_rel0] = lin_MCK_consistent(tY{r}, aY{r}, M, Cstr, K);
        w = (tY{r} >= t5y(r) & tY{r} <= t95y(r));
        drift0 = x0(:,2:end) - x0(:,1:end-1);
        idr0   = abs(drift0(w,:)) ./ (ones(sum(w),1) * h_story(:)');
        refIDR.Y.d(r) = max(idr0,[],'all');
        
        a_abs0 = a_rel0 + aY{r}(:) * ones(1, size(a_rel0,2));
        A_i = zeros(1, size(a_abs0,2));
        for i=1:size(a_abs0,2)
            ai = a_abs0(w,i);
            rmsv = sqrt(mean(ai.^2,'omitnan'));
            p95  = prctile(abs(ai),95);
            A_i(i) = rmsv + w_p * p95;
        end
        refACC.Y.a(r) = max(A_i);
    end
end

[J1, ~] = compute_J1_IDR_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, design_set, x_ga, refIDR);

[J2, ~] = compute_J2_ACC_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, env, design_set, x_ga, refACC);
end

%% ---------------------------------------------------------------------
function key = hash_psa_prep(t_cells, a_cells, params)
%HASH_PSA_PREP Create hash for PSA preparation inputs
buf = [];
for i = 1:numel(t_cells)
    if ~isempty(t_cells{i}), buf = [buf; t_cells{i}(:)]; end
    if ~isempty(a_cells{i}), buf = [buf; a_cells{i}(:)]; end
end
buf = [buf; params(:)];
key = md5(num2str(buf(:)', '%.16g,'));
end










