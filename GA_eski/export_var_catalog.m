function export_var_catalog()
%EXPORT_VAR_CATALOG Collects key model variables and writes out/var_catalog.xlsx
% Categories: 'orifis', 'termal', 'hat_atalet', 'arias', 'cms', 'diger'
% Run: export_var_catalog;  -> writes out/var_catalog.xlsx

    if ~exist('out','dir'), mkdir('out'); end
    xfile = fullfile('out','var_catalog.xlsx');

    E = {};

    % ---- Orifis (orifice / flow) ----
    add('orifis','geom.d_o','Orifis Ã§apÄ± [m]','viscous.m','geom.d_o');
    add('orifis','geom.Lori','Orifis uzunluÄŸu [m]','viscous.m','geom.Lori');
    add('orifis','orf.n_orf','Orifis adedi [-]','viscous.m','orf.n_orf');
    add('orifis','orf.Cd0','DÃ¼ÅŸÃ¼k Re sÃ¼rÃ¼kleme katsayÄ±sÄ± [-]','viscous.m','orf.Cd0');
    add('orifis','orf.CdInf','YÃ¼ksek Re sÃ¼rÃ¼kleme katsayÄ±sÄ± [-]','viscous.m','orf.CdInf');
    add('orifis','orf.Rec','Reynolds geÃ§iÅŸ eÅŸiÄŸi [-]','viscous.m','orf.Rec');
    add('orifis','orf.p_exp','Cd(Re) eÄŸri Ã¼ssÃ¼ [-]','viscous.m','orf.p_exp');
    add('orifis','orf.cav_sf','Kavitasyon gÃ¼venlik Ã§arpanÄ± [-]','viscous.m','orf.cav_sf');
    add('orifis','num.dP_cap','BasÄ±nÃ§ sÄ±nÄ±rlandÄ±rma (cap) [Pa]','viscous.m','num.dP_cap');
    add('orifis','cfg.on.CdRe','Cd(Re) modeli aÃ§Ä±k/kapa','viscous.m','cfg.on.CdRe');
    add('orifis','cfg.on.Rkv','Kareli kayÄ±p (kv) aÃ§Ä±k/kapa','viscous.m','cfg.on.Rkv');
    add('orifis','cfg.on.Rlam','Laminer kayÄ±p (R_lam) aÃ§Ä±k/kapa','viscous.m','cfg.on.Rlam');
    add('orifis','cfg.on.Qsat','AkÄ±ÅŸ satÃ¼rasyonu (tanh) aÃ§Ä±k/kapa','viscous.m','cfg.on.Qsat');
    add('orifis','num.Qcap_big','Referans bÃ¼yÃ¼k Q kapasitesi [m^3/s] (tÃ¼retilmiÅŸ)','viscous.m','num.Qcap_big');

    % ---- Termal ----
    add('termal','therm.mu_ref','Referans viskozite [PaÂ·s]','viscous.m','therm.mu_ref');
    add('termal','therm.b_mu','Î¼(T) eÄŸimi [1/Â°C]','viscous.m','therm.b_mu');
    add('termal','therm.beta0','Referans hacim modÃ¼lÃ¼ [Pa]','viscous.m','therm.beta0');
    add('termal','therm.b_beta','Î²(T) eÄŸimi [1/Â°C]','viscous.m','therm.b_beta');
    add('termal','therm.hA_os','YaÄŸ-Ã§elik Ä±sÄ± geÃ§iÅŸ katsayÄ±sÄ± [W/K]','viscous.m','therm.hA_os');
    add('termal','therm.hA_o_env','YaÄŸ-Ã§evre Ä±sÄ± geÃ§iÅŸi [W/K]','mck_with_damper_adv.m','therm.hA_o_env');
    add('termal','therm.hA_s_env','Ã‡elik-Ã§evre Ä±sÄ± geÃ§iÅŸi [W/K]','mck_with_damper_adv.m','therm.hA_s_env');
    add('termal','therm.resFactor','YaÄŸ hacim Ã¶lÃ§eÄŸi (V_oil) Ã§arpanÄ± [-]','viscous.m','therm.resFactor');
    add('termal','therm.cp_oil','YaÄŸÄ±n c_p [J/(kgÂ·K)]','viscous.m','therm.cp_oil');
    add('termal','therm.cp_steel','Ã‡eliÄŸin c_p [J/(kgÂ·K)]','viscous.m','therm.cp_steel');
    add('termal','therm.rho_ref','Referans yoÄŸunluk [kg/m^3]','viscous.m','therm.rho_ref');
    add('termal','therm.alpha_rho','YoÄŸunluk sÄ±caklÄ±k katsayÄ±sÄ± [1/Â°C]','viscous.m','therm.alpha_rho');
    add('termal','therm.T_ref_C','Referans sÄ±caklÄ±k [Â°C]','mck_with_damper_adv.m','therm.T_ref_C');
    add('termal','therm.T_env_C','Ã‡evre sÄ±caklÄ±ÄŸÄ± [Â°C]','mck_with_damper_adv.m','therm.T_env_C');
    add('termal','therm.T0_C','BaÅŸlangÄ±Ã§ yaÄŸ sÄ±caklÄ±ÄŸÄ± [Â°C]','mck_with_damper_adv.m','therm.T0_C');
    add('termal','therm.Ts0_C','BaÅŸlangÄ±Ã§ Ã§elik sÄ±caklÄ±ÄŸÄ± [Â°C]','mck_with_damper_adv.m','therm.Ts0_C');
    add('termal','cfg.use_thermal','Termal model aÃ§Ä±k/kapa','viscous.m','cfg.use_thermal');
    add('termal','cfg.on.mu_floor','Î¼ alt sÄ±nÄ±r korumasÄ± aÃ§Ä±k/kapa','viscous.m','cfg.on.mu_floor');

    % ---- Hat atalet / hidrolik Ã§evre ----
    add('hat_atalet','hyd.Lh','Hat atalet parametresi [kg/m^4 ~ PaÂ·s^2/m^3]','viscous.m','hyd.Lh');
    add('hat_atalet','hyd.K_leak','SÄ±zÄ±ntÄ± iletkenliÄŸi [m^3/(sÂ·Pa)]','viscous.m','hyd.K_leak');
    add('hat_atalet','hyd.Vmin_fac','Minimum hacim skala faktÃ¶rÃ¼ [-]','viscous.m','hyd.Vmin_fac');
    add('hat_atalet','hyd.n_parallel','Paralel damper sayÄ±sÄ± [-]','mck_with_damper_adv.m','hyd.n_parallel');
    add('hat_atalet','hyd.V0','Hazne referans hacmi (tÃ¼retilmiÅŸ) [m^3]','viscous.m','hyd.V0');
    add('hat_atalet','cfg.on.hyd_inertia','Hat atalet modeli aÃ§Ä±k/kapa','viscous.m','cfg.on.hyd_inertia');
    add('hat_atalet','cfg.on.leak','SÄ±zÄ±ntÄ± modeli aÃ§Ä±k/kapa','viscous.m','cfg.on.leak');
    add('hat_atalet','cfg.on.cavitation','Kavitasyon modeli aÃ§Ä±k/kapa','viscous.m','cfg.on.cavitation');
    add('hat_atalet','num.softmin_eps','YumuÅŸak minimum iÃ§in epsilon [Pa]','mck_with_damper_adv.m','num.softmin_eps');

    % ---- Arias penceresi ----
    add('arias','pp.on.arias','Arias penceresi etkinliÄŸi','viscous.m','pp.on.arias');
    add('arias','pp.tail_sec','SimÃ¼lasyon kuyruÄŸu sÃ¼resi [s]','viscous.m','pp.tail_sec');
    add('arias','obj.use_arias_window','Hedef metrikte Arias penceresi kullan','viscous.m','obj.use_arias_window');
    add('arias','obj.window_source','Pencere kaynaÄŸÄ±: raw/scaled/same','viscous.m','obj.window_source');

    % ---- CMS / PSA / yoÄŸunluk eÅŸitleme ----
    add('cms','pp.on.intensity','PSA band ort. ÅŸiddet normalizasyonu','viscous.m','pp.on.intensity');
    add('cms','pp.on.CMS','CMS hedefi aÃ§Ä±k/kapa','viscous.m','pp.on.CMS');
    add('cms','pp.gammaCMS','Band vs CMS aÄŸÄ±rlÄ±ÄŸÄ± [0..1]','viscous.m','pp.gammaCMS');
    add('cms','pp.PSA.zeta','PSA tek serbestlik sÃ¶nÃ¼m oranÄ± [-]','viscous.m','pp.PSA.zeta');
    add('cms','pp.PSA.band_fac','PSA T1 band faktÃ¶rÃ¼ [min max]','viscous.m','pp.PSA.band_fac');
    add('cms','pp.PSA.Np_band','Bant iÃ§i periyot sayÄ±sÄ± [-]','viscous.m','pp.PSA.Np_band');
    add('cms','pp.PSA.downsample_enabled','PSA iÃ§in downsample aÃ§/kapa','viscous.m','pp.PSA.downsample_enabled');
    add('cms','pp.PSA.downsample_dt','Downsample hedef dt [s]','viscous.m','pp.PSA.downsample_dt');
    add('cms','pp.PSA.use_parfor','PSA paralel hesap aÃ§/kapa','viscous.m','pp.PSA.use_parfor');

    % ---- DiÄŸer (geometri, yay, PF, amaÃ§) ----
    add('diger','geom.Dp','Piston Ã§apÄ± [m]','viscous.m','geom.Dp');
    add('diger','geom.Lgap','BoÅŸluk (gap) [m]','viscous.m','geom.Lgap');
    add('diger','geom.Kd','Hidrolik eÅŸdeÄŸer rijitlik Ã¶lÃ§eÄŸi [Pa]','viscous.m','geom.Kd');
    add('diger','geom.Ebody','GÃ¶vde elastisite modÃ¼lÃ¼ [Pa]','viscous.m','geom.Ebody');
    add('diger','sh.d_w','Yay tel Ã§apÄ± [m]','viscous.m','sh.d_w');
    add('diger','sh.D_m','Yay ortalama Ã§apÄ± [m]','viscous.m','sh.D_m');
    add('diger','sh.n_turn','Yay sarÄ±m sayÄ±sÄ± [-]','viscous.m','sh.n_turn');
    add('diger','sh.G','Yay kesme modÃ¼lÃ¼ [Pa]','simulate.m','sh.G');
    add('diger','cfg.PF.mode','BasÄ±nÃ§-kuvvet mod: ramp','viscous.m','cfg.PF.mode');
    add('diger','cfg.PF.auto_t_on','PF t_on otomatik ayarla','viscous.m','cfg.PF.auto_t_on');
    add('diger','cfg.PF.t_on','PF baÅŸlama zamanÄ± [s]','viscous.m','cfg.PF.t_on');
    add('diger','cfg.PF.tau','PF zaman sabiti [s]','viscous.m','cfg.PF.tau');
    add('diger','cfg.PF.k','PF softness param [s]','pf_weight.m','cfg.PF.k');
    add('diger','cfg.PF.gain','PF kazanÃ§ [-]','viscous.m','cfg.PF.gain');
    add('diger','cfg.on.pf_resistive_only','Sadece rezistif PF bileÅŸeni','viscous.m','cfg.on.pf_resistive_only');
    add('diger','cfg.on.pressure_ode','BasÄ±nÃ§ ODE aÃ§Ä±k/kapa','viscous.m','cfg.on.pressure_ode');
    add('diger','cfg.on.pressure_force','BasÄ±nÃ§ kuvvet katkÄ±sÄ± aÃ§Ä±k/kapa','viscous.m','cfg.on.pressure_force');
    add('diger','obj.idx_disp_story','Ä°zlenen kat (yer deÄŸiÅŸtirme)','viscous.m','obj.idx_disp_story');
    add('diger','obj.idx_acc_story','Ä°zlenen kat (ivme)','viscous.m','obj.idx_acc_story');
    add('diger','obj.acc_metric','Ä°vme metrik seÃ§imi','viscous.m','obj.acc_metric');
    add('diger','obj.p95_penalty_w','Ä°vme p95 aÄŸÄ±rlÄ±k Ã§arpanÄ±','viscous.m','obj.p95_penalty_w');

    % ---- Yaz ----
    cats_all = cellfun(@(s) char(s.category), E, 'UniformOutput', false);
    cats = unique(cats_all); % cell array of char vectors
    for ci = 1:numel(cats)
        cat = cats{ci};
        rows = E(strcmp(cat, cats_all));
        T = struct2table(cell2mat(rows(:)), 'AsArray', true);
        % Kolon sÄ±rasÄ±
        T = T(:, {'category','name','path','description','source'});
        writetable(T, xfile, 'Sheet', cat, 'WriteMode','overwritesheet');
    end
    fprintf('YazÄ±ldÄ±: %s (sayfalar: %s)\n', xfile, strjoin(cats, ', '));

    function add(category, name, description, source, path)
        S = struct('category',string(category), 'name',string(name), ...
                   'description',string(description), 'source',string(source), ...
                   'path',string(path));
        E{end+1} = S; %#ok<AGROW>
    end
end

