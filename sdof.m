%% SDOF analysis using GA best parameters
% Runs a single simulation without optimization using the best design
% parameters (xbest) stored in ga_front.csv. Produces displacement and
% acceleration histories for the 10th floor along with peak metrics.

clear; clc; close all;

%% Load base parameters (embedded defaults)
params = default_params();

M = params.M;        C0 = params.C0;       K = params.K;
n = params.n;        story_height = params.story_height;
T1 = params.T1;

Dp = params.Dp;      Lgap = params.Lgap;   d_w = params.d_w;
D_m = params.D_m;    n_turn = params.n_turn;
mu_ref = params.mu_ref;  Lori = params.Lori;
n_orf = params.n_orf;    rho = params.rho;
orf = params.orf;        cfg = params.cfg;
thermal = params.thermal;
T0_C = params.T0_C;  T_ref_C = params.T_ref_C;  b_mu = params.b_mu;
cp_oil = params.cp_oil;   cp_steel = params.cp_steel;
steel_to_oil_mass_ratio = params.steel_to_oil_mass_ratio;
n_dampers_per_story = params.n_dampers_per_story;
story_mask = params.story_mask;
k_hyd_single = params.k_hyd_single;
k_p_single = params.k_p;
Kd = params.Kd; Ebody = params.Ebody; Gsh = params.Gsh;

Ap = params.Ap;   Ao = params.Ao;   k_sd = params.k_sd;
c_lam0 = params.c_lam0;  c_lam_min_frac = params.c_lam_min_frac;
c_lam_min = params.c_lam_min; c_lam_max = params.c_lam_max;
Qcap_big = params.Qcap_big; dP_max = params.dP_max;
V_oil_per = params.V_oil_per;

%% Override with GA best parameters from GA output
try
    tbl = readtable('ga_front.csv');
    xb = tbl(1,:);

    d_o = xb.d_o_mm/1000;      % [m]
    n_orf = xb.n_orf;
    cfg.PF.tau  = xb.PF_tau;
    cfg.PF.gain = xb.PF_gain;
    cfg.PF.t_on = xb.PF_t_on;
    orf.Cd0  = xb.Cd0;
    orf.CdInf= xb.CdInf;
    orf.p_exp= xb.p_exp;
    Lori  = xb.Lori_mm/1000;    % [m]
    thermal.hA_W_perK = xb.hA_W_perK;
    Dp    = xb.Dp_mm/1000;      % [m]
    d_w   = xb.d_w_mm/1000;     % [m]
    D_m   = xb.D_m_mm/1000;     % [m]
    n_turn= xb.n_turn;
    mu_ref= xb.mu_ref;

    orf.d_o = d_o;

    params_override = struct('Dp',Dp,'d_w',d_w,'D_m',D_m,'n_turn',n_turn, ...
        'mu_ref',mu_ref,'Lori',Lori,'Lgap',Lgap,'Kd',Kd,'Ebody',Ebody, ...
        'Gsh',Gsh,'orf',orf,'n_orf',n_orf,'rho',rho, ...
        'c_lam_min_frac',c_lam_min_frac);

    params_override = build_params_local(params_override);

    Ap = params_override.Ap;
    Ao = params_override.Ao;
    k_sd = params_override.k_sd;
    k_hyd_single = params_override.k_hyd_single;
    k_p_single = params_override.k_p;
    c_lam0 = params_override.c_lam0;
    c_lam_min = params_override.c_lam_min;
    Qcap_big = params_override.Qcap_big;
    dP_max = params_override.dP_max;
    V_oil_per = params_override.V_oil_per;
    orf = params_override.orf;
catch ME
    warning('Unable to apply GA parameters: %s', ME.message);
end

%% Ground motion input (scaled by T1)
[~, recs] = load_ground_motions_local(T1);
rec = recs(1);  t = rec.t; ag = rec.ag;
win = make_arias_window_local(t, ag); t5 = win.t5; t95 = win.t95;

%% Damperless and damper responses
[x0,a_rel0] = solve_linear_mck_local(t, ag, M, C0, K);
[x_d,a_d,diag_d] = mck_with_damper_local(t, ag, M, C0, K, k_sd, c_lam0, Lori, ...
    orf, rho, Ap, Ao, k_hyd_single, k_p_single, mu_ref, thermal, T0_C, T_ref_C, b_mu, ...
    c_lam_min, c_lam_max, Lgap, cp_oil, cp_steel, steel_to_oil_mass_ratio, ...
    story_mask, n_dampers_per_story, V_oil_per, cfg, dP_max);

x10_0 = x0(:,10);
x10_d = x_d(:,10);
a10_0 = a_rel0(:,10) + ag;
a10_d = a_d(:,10) + ag;

%% Assemble equivalent damping/stiffness matrices for modal check
nStories = n - 1;
mask = story_mask(:); if numel(mask)==1, mask = mask*ones(nStories,1); end
ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
multi = (mask .* ndps);
Kadd = zeros(n); C_add = zeros(n);
c_diag = diag_d.c_lam_inst(end,:);
if numel(c_diag)==1
    c_diag = repmat(c_diag, 1, nStories);
end
for i = 1:nStories
    idx = [i,i+1];
    k_eq = k_sd * multi(i);
    c_eq = c_diag(i);
    kM  = k_eq * [1 -1; -1 1];
    cM  = c_eq * [1 -1; -1 1];
    Kadd(idx,idx) = Kadd(idx,idx) + kM;
    C_add(idx,idx)= C_add(idx,idx) + cM;
end
K_tot = K + Kadd;
C_d = C0 + C_add;
[V,D] = eig(K_tot,M); [w2,ord] = sort(diag(D),'ascend');
phi1 = V(:,ord(1)); w1 = sqrt(w2(1));
normM = phi1.' * M * phi1;

%% Plots: 10th floor displacement and acceleration, plus IDR
figure('Name','10. Kat yer değiştirme — ham ivme (ODE-only)','Color','w');
plot(t, x10_0,'k','LineWidth',1.4); hold on;
plot(t, x10_d,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('x10(t) [m]');
title(sprintf('10-Kat | T1=%.3f s | Arias [%.3f, %.3f] s', T1, t5, t95));
legend('Dampersiz','Damperli','Location','best');

figure('Name','10. Kat mutlak ivme','Color','w');
plot(t, a10_0,'k','LineWidth',1.4); hold on;
plot(t, a10_d,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('a10abs(t) [m/s^2]');
legend('Dampersiz','Damperli','Location','best');

drift0 = x0(:,2:end) - x0(:,1:end-1);
drift_d = x_d(:,2:end) - x_d(:,1:end-1);
idx = win.idx;
IDR0 = max(abs(drift0(idx,:)))./story_height;
IDR_d = max(abs(drift_d(idx,:)))./story_height;
story_ids = 1:(n-1);
figure('Name','Maksimum IDR','Color','w');
plot(story_ids, IDR0,'k-o','LineWidth',1.4); hold on;
plot(story_ids, IDR_d,'r-d','LineWidth',1.0);
grid on; xlabel('Kat'); ylabel('Maks IDR [Delta x/h]');
legend('Dampersiz','Damperli','Location','best');

%% Summary printout
zeta0 = (phi1.' * C0 * phi1) / (2*w1*normM);
zeta_d = (phi1.' * C_d * phi1) / (2*w1*normM);

idx = win.idx;
x10_max_0 = max(abs(x10_0(idx)));
x10_max_damperli = max(abs(x10_d(idx)));
a10abs_max_0 = max(abs(a10_0(idx)));
a10abs_max_damperli = max(abs(a10_d(idx)));
IDR_max_0 = max(IDR0);
IDR_max_damperli = max(IDR_d);

fprintf(['Self-check zeta1: %.3f%% (dampersiz) vs %.3f%% (damperli); ' ...
         'x10_{max}0=%.4g m; x10_{max}d=%.4g m; ' ...
         'a10abs_{max}0=%.4g m/s^2; a10abs_{max}d=%.4g m/s^2; ' ...
         'IDR_{max}0=%.4g; IDR_{max}d=%.4g\n'], ...
        100*zeta0, 100*zeta_d, x10_max_0, x10_max_damperli, ...
        a10abs_max_0, a10abs_max_damperli, IDR_max_0, IDR_max_damperli);

%% Local helpers
function params = default_params()
    % Inline version of parametreler.m providing baseline structural and
    % damper parameters together with derived quantities.
    n  = 10;
    m  = 360e3 * ones(n,1);
    k  = 6.5e8 * ones(n,1);
    c  = 6.2e6 * ones(n,1);
    story_height = 3.0;

    M  = diag(m);
    K  = zeros(n); C0 = zeros(n);
    for ii = 1:n
        kL = k(ii); cL = c(ii);
        if ii < n
            kU = k(ii+1); cU = c(ii+1);
        else
            kU = 0; cU = 0;
        end
        K(ii,ii) = kL + kU; C0(ii,ii) = cL + cU;
        if ii > 1
            K(ii,ii-1) = -kL; C0(ii,ii-1) = -cL;
        end
        if ii < n
            K(ii,ii+1) = -kU; C0(ii,ii+1) = -cU;
        end
    end

    [~,D] = eig(K,M);
    w  = sqrt(sort(diag(D),'ascend'));
    T1 = 2*pi/w(1);

    params = struct();
    params.n = n;
    params.M = M; params.C0 = C0; params.K = K;
    params.story_height = story_height;
    params.T1 = T1;

    params.Dp   = 0.125;
    params.Lgap = 0.055;
    % Orifis + hazne parametreleri: smooth_eps0 yumuşatma, reservoir ise
    % lineer yay/damper katsayıları ile minimum basıncı temsil eder.
    % p_cav_eff buhar basıncı tabanını belirler; varsayılan olarak p_amb'e
    % (veya orf.p_supply_cap ile belirtilen başka bir değere) kadar olan
    % aralıkta çekiş kuvvetlerini sınırlar ve hazne basıncı ambientteyken
    % dahi sınırlı bir emiş düşümü sağlar.
    params.orf = struct('Cd0',0.61,'CdInf',0.80,'Rec',3000,'p_exp',1.1, ...
        'p_amb',1.0e5,'p_cav_eff',2.0e3,'cav_sf',0.90,'d_o',3.0e-3,'veps',0.10, ...
        'smooth_eps0',5.0e2,'reservoir',struct('p_min',2.0e3,'k_lin',4.0e8, ...
        'c_lin',2.0e5,'V_ref',2.5e-4));
    params.Lori = 0.10;
    params.mu_ref = 0.9;

    params.Kd = 1.5e9;   % Yağın hacimsel modülü [Pa]
    params.Ebody = 2.1e11;
    params.Gsh = 79e9;
    params.d_w = 12e-3;
    params.D_m = 80e-3;
    params.n_turn = 8;

    params.rho = 850;
    params.n_orf = 6;

    params.T0_C = 25;
    params.T_ref_C = 25;
    params.b_mu = -0.013;
    params.thermal = struct('hA_W_perK',450,'T_env_C',25,'max_iter',3, ...
        'tol_K',0.5,'relax',0.5,'T_floor_C',-40,'T_ceil_C',180);

    params.steel_to_oil_mass_ratio = 1.5;
    params.n_dampers_per_story = 1;
    params.story_mask = ones(n-1,1);
    params.cp_oil = 1800;
    params.cp_steel = 500;

    params.c_lam_min_frac = 0.05;

    cfg = struct();
    % cfg.PF akümülatör/valf fiziğini tarif eder; pf_weight_local artık yoktur.
    cfg.PF = struct('mode','accumulator', ...
        'precharge_Pa',1.05e5,'gas_gamma',1.2,'V_gas_ref',2.5e-4, ...
        'bleed_coeff',2.0e5,'notes','Fiziksel akümülatör modeli; kuvvet ölçekleme yok');
    cfg.on = struct('pressure_force',true,'mu_floor',false);
    cfg.compat_simple = false;
    cfg.num = struct('mu_min_phys',0.6,'dP_cap',NaN);
    params.cfg = cfg;

    params = build_params_local(params);
end

function params = build_params_local(params)
    params = recompute_damper_params_local(params);
    params.Qcap_big = 0;
    params.dP_max = 0;
    if isfield(params,'orf') && isfield(params.orf,'CdInf') && ...
            isfield(params,'Ao') && isfield(params,'rho') && isfield(params,'Ap')
        p_cav_floor = max(params.orf.p_cav_eff, 0);
        p_supply_cap = getfield_default_local(params.orf,'p_supply_cap', params.orf.p_amb);
        p_supply_cap = max(p_supply_cap, p_cav_floor);
        dP_cav = max(p_supply_cap - p_cav_floor, 0) * params.orf.cav_sf;
        stroke_ref = getfield_default_local(params,'stroke_ref',params.Lgap);
        V_ref = max(getfield_default_local(params,'V_oil_per',params.Ap*(2*params.Lgap)), 1e-12);
        dV_struct = max(params.Ap * stroke_ref, 0);
        dP_struct = params.Kd * dV_struct / V_ref;
        dP_max = min(dP_struct, dP_cav);
        cfg_dP = NaN;
        if isfield(params,'cfg') && isfield(params.cfg,'num')
            cfg_dP = getfield_default_local(params.cfg.num,'dP_cap',NaN);
        end
        if isfinite(cfg_dP)
            dP_max = min(dP_max, cfg_dP);
        end
        params.dP_max = max(dP_max, 0);
        if params.dP_max > 0
            params.Qcap_big = max(params.orf.CdInf * params.Ao, 1e-12) * ...
                sqrt(2*params.dP_max / params.rho);
        end
    end

    if isfield(params,'c_lam0') && isfield(params,'orf') && ...
            isfield(params.orf,'d_o') && isfield(params,'Lori') && ...
            isfield(params,'Ap')
        mu_ref = getfield_default_local(params,'mu_ref',1);
        mu_min_ratio = NaN;
        mu_max_ratio = NaN;
        if isfield(params,'cfg') && isfield(params.cfg,'num')
            mu_min_ratio = getfield_default_local(params.cfg.num,'mu_min_phys',NaN);
            mu_max_ratio = getfield_default_local(params.cfg.num,'mu_max_phys',NaN);
        end
        nd = getfield_default_local(params,'n_parallel',1);
        mu_min = mu_ref;
        if isfinite(mu_min_ratio)
            mu_min = max(mu_ref * mu_min_ratio, 1e-6);
        end
        mu_max = mu_ref;
        if isfinite(mu_max_ratio)
            mu_max = max(mu_ref * mu_max_ratio, mu_min);
        end
        c_from_mu = @(mu) nd * (128/pi) * mu * params.Lori * (params.Ap.^2) / (params.orf.d_o.^4);
        c_lam_min_phys = c_from_mu(mu_min);
        c_lam_max_phys = c_from_mu(mu_max);
        frac_floor = getfield_default_local(params,'c_lam_min_frac',0);
        params.c_lam_min = max(c_lam_min_phys, frac_floor * params.c_lam0);
        params.c_lam_max = max(c_lam_max_phys, params.c_lam0);
    end
end

function params = recompute_damper_params_local(params)
    if ~isstruct(params), return; end

    if isfield(params,'Dp_mm'),    params.Dp   = params.Dp_mm/1000; end
    if isfield(params,'d_w_mm'),   params.d_w  = params.d_w_mm/1000; end
    if isfield(params,'D_m_mm'),   params.D_m  = params.D_m_mm/1000; end
    if isfield(params,'Lori_mm'),  params.Lori = params.Lori_mm/1000; end
    if isfield(params,'orf') && isfield(params.orf,'d_o_mm')
        params.orf.d_o = params.orf.d_o_mm/1000;
    end

    req = {'Dp','d_w','D_m','n_turn','mu_ref','Lori','Lgap','Kd','Ebody','Gsh'};
    if ~all(isfield(params,req)) || ~isfield(params,'orf') || ~isfield(params.orf,'d_o')
        return;
    end

    if ~isfield(params,'n_orf') && isfield(params,'orf') && isfield(params.orf,'n_orf')
        params.n_orf = params.orf.n_orf;
    end
    if ~isfield(params,'n_orf'), params.n_orf = 1; end

    nd = 1;
    if isfield(params,'nd') && isfinite(params.nd)
        nd = params.nd;
    elseif isfield(params,'n_parallel') && isfinite(params.n_parallel)
        nd = params.n_parallel;
    elseif isfield(params,'n_dampers_per_story')
        nds = params.n_dampers_per_story;
        if isnumeric(nds)
            if isscalar(nds), nd = nds; else, nd = max(1, round(max(nds(:)))); end
        end
    end
    nd = max(1, round(nd));

    Ap_single = pi * params.Dp^2 / 4;
    Ao_single = pi * params.orf.d_o^2 / 4;
    Ao_damper = params.n_orf * Ao_single;
    Ap_eff = nd * Ap_single;
    Ao_eff = nd * Ao_damper;

    V_oil_geom_single = Ap_single*(2*params.Lgap) + params.n_orf*Ao_single*params.Lori;
    V_oil_extra = getfield_default_local(params,'V_oil_extra',0);
    V_oil_single = max(V_oil_geom_single + V_oil_extra, 0);

    rho_loc = getfield_default_local(params,'rho',850);
    Lh = rho_loc * params.Lori / max(Ao_eff, 1e-12);

    k_bulk_single = params.Kd * Ap_single^2 / max(V_oil_single, 1e-12);
    k_body_single = params.Ebody * Ap_single / max(params.Lgap, 1e-12);
    inv_k_bulk = inf;
    if isfinite(k_bulk_single) && k_bulk_single > 0
        inv_k_bulk = 1 / k_bulk_single;
    end
    inv_k_body = inf;
    if isfinite(k_body_single) && k_body_single > 0
        inv_k_body = 1 / k_body_single;
    end
    inv_sum = inv_k_bulk + inv_k_body;
    if inv_sum <= 0
        k_hyd_single = 0;
    else
        k_hyd_single = 1 / inv_sum;
    end
    k_p = params.Gsh * params.d_w^4 / (8 * params.n_turn * params.D_m^3);
    k_sd_simple = k_hyd_single + k_p;
    k_sd_adv    = nd * (k_hyd_single + k_p);

    c_lam_single = (128/pi) * params.mu_ref * params.Lori * Ap_single^2 / (params.orf.d_o^4);
    c_lam0 = nd * c_lam_single;

    params.V_oil_per = V_oil_single;

    params.Ap = Ap_single;
    params.Ao = Ao_damper;
    params.Ap_eff = Ap_eff;
    params.Ao_eff = Ao_eff;
    params.Lh = Lh;
    params.V_oil_single = V_oil_single;
    params.V_oil_total = V_oil_single * nd;
    params.n_parallel = nd;
    params.k_p = k_p;
    params.k_hyd_single = k_hyd_single;
    params.k_bulk_single = k_bulk_single;
    params.k_sd_simple = k_sd_simple;
    params.k_sd_adv = k_sd_adv;
    params.k_sd = k_sd_adv;
    params.c_lam0 = c_lam0;
    if isfinite(k_hyd_single) && k_hyd_single > 0
        params.hydraulic_compliance_single = 1 ./ k_hyd_single;
    else
        params.hydraulic_compliance_single = inf;
    end
end

function val = getfield_default_local(S, fname, defaultVal)
    if ~isstruct(S) || ~isfield(S,fname) || isempty(S.(fname))
        val = defaultVal; return;
    end
    valRaw = S.(fname);
    if isnumeric(valRaw)
        if isempty(valRaw) || any(~isfinite(valRaw(:)))
            val = defaultVal;
        else
            val = valRaw;
        end
    else
        val = valRaw;
    end
end

function y = softmin_local(a,b,epsm)
    if nargin < 3 || isempty(epsm)
        epsm = 1e5;
    end
    y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
end

function win = make_arias_window_local(t, ag, varargin)
    p = inputParser;
    p.addParameter('p1',0.05,@(x)isscalar(x) && x>=0 && x<=1);
    p.addParameter('p2',0.95,@(x)isscalar(x) && x>=0 && x<=1);
    p.addParameter('pad',0.5,@(x)isscalar(x) && x>=0);
    p.parse(varargin{:});
    p1 = p.Results.p1; p2 = p.Results.p2; pad = p.Results.pad;

    IA = cumtrapz(t, ag.^2);
    IA_tot = IA(end);
    if ~(isfinite(IA_tot)) || IA_tot <= eps
        t_start = t(1); t_end = t(end);
        idx = true(size(t));
        win = struct('t5',t_start,'t95',t_end,'pad',0, ...
                     't_start',t_start,'t_end',t_end,'idx',idx, ...
                     'coverage',1.0,'flag_low_arias',true);
        return;
    end

    IA_norm = IA / IA_tot;
    [t_unique, iu] = unique(IA_norm);
    interp_t = t(iu);
    t5  = interp1(t_unique, interp_t, p1, 'linear');
    t95 = interp1(t_unique, interp_t, p2, 'linear');
    dur = t95 - t5;
    if dur < 5, pad = 0.25; end
    t_start = max(t(1),  t5  - pad);
    t_end   = min(t(end), t95 + pad);
    idx = (t >= t_start) & (t <= t_end);
    coverage = trapz(t(idx), ag(idx).^2) / IA_tot;
    flag_low_arias = coverage < 0.90;
    win = struct('t5',t5,'t95',t95,'pad',pad, ...
                 't_start',t_start,'t_end',t_end,'idx',idx, ...
                 'coverage',coverage,'flag_low_arias',flag_low_arias);
end

function [records, scaled, meta] = load_ground_motions_local(T1, opts)
    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'verbose'), opts.verbose = true; end
    hp_cut   = getfield_default_local(opts,'hp_cut',0.05);
    IM_mode  = getfield_default_local(opts,'IM_mode','band');
    band_fac = getfield_default_local(opts,'band_fac',[0.8 1.2]);
    band_N   = getfield_default_local(opts,'band_N',21);
    s_bounds = getfield_default_local(opts,'s_bounds',[0.2 2.2]);

    if exist('acc_matrix.mat','file')
        raw = load('acc_matrix.mat');
    else
        warning('load_ground_motions_local:synthetic', ...
            'acc_matrix.mat bulunamadı; sentetik tek kayıt kullanılacak.');
        t_syn = (0:0.01:40)';
        ag_syn = 0.35 * sin(2*pi*1.2*t_syn) .* exp(-0.05*t_syn);
        raw = struct('synthetic',[t_syn ag_syn]);
    end
    fn  = fieldnames(raw);
    records = struct('name',{},'t',{},'ag',{},'dt',{},'duration',{},'PGA',{},'PGV',{},'IM',{},'scale',{});

    for kk = 1:numel(fn)
        A = raw.(fn{kk});
        t  = A(:,1);   ag = A(:,2);
        [t,iu] = unique(t,'stable'); ag = ag(iu);
        dt = median(diff(t));
        assert(max(abs(diff(t) - dt)) < 1e-6, 'Zaman örnekleme aralığı düzensiz');
        t  = (t(1):dt:t(end)).';
        ag = interp1(A(:,1), A(:,2), t, 'linear');
        assert(max(abs(ag)) < 100, 'İvme büyüklüğü birim hatasına işaret ediyor');
        ag = detrend(ag,0);
        ag = detrend(ag,1);
        try
            Fs = 1/dt;
            Wn = hp_cut/(Fs/2);
            [b,a] = butter(2,Wn,'high');
            ag = filtfilt(b,a,ag);
        catch
        end
        duration = t(end) - t(1);
        v  = cumtrapz(t,ag);
        PGA = max(abs(ag));
        PGV = max(abs(v));
        records(end+1) = struct('name',fn{kk},'t',t,'ag',ag,'dt',dt, ...
                                 'duration',duration,'PGA',PGA,'PGV',PGV, ...
                                 'IM',[],'scale',1); %#ok<AGROW>
    end

    if getfield_default_local(opts,'verbose',true)
        fprintf('Toplam %d zemin hareketi kaydı yüklendi:\n', numel(records));
        for kk = 1:numel(records)
            r = records(kk);
            fprintf('%2d) %-12s dt=%6.4f s dur=%6.2f s PGA=%7.3f PGV=%7.3f\n', ...
                kk, r.name, r.dt, r.duration, r.PGA, r.PGV);
        end
    end

    scaled = [];
    meta = struct();

    if nargin >= 1 && ~isempty(T1)
        for kk = 1:numel(records)
            records(kk).IM = compute_IM_local(records(kk).t, records(kk).ag, IM_mode, T1, band_fac, band_N);
        end

        IM = [records.IM];
        IM_low  = max(s_bounds(1)*IM);
        IM_high = min(s_bounds(2)*IM);
        targetIM0 = median(IM);
        targetIM  = min(max(targetIM0, IM_low), IM_high);
        doClip = (IM_low > IM_high);

        max_iter = 3; dropped = {};
        for it = 1:max_iter
            IM_low  = max(s_bounds(1)*IM);
            IM_high = min(s_bounds(2)*IM);
            if IM_low <= IM_high, break; end
            [~,idx] = max(abs(log(IM) - median(log(IM))));
            dropped{end+1} = records(idx).name; %#ok<AGROW>
            records(idx) = [];
            IM(idx) = [];
        end
        if ~isempty(dropped) && getfield_default_local(opts,'verbose',true)
            fprintf('TRIM: ayıklanan uç değerler = %s\n', strjoin(dropped,', '));
        end
        IM_low  = max(s_bounds(1)*IM);
        IM_high = min(s_bounds(2)*IM);
        targetIM0 = median(IM);
        targetIM  = min(max(targetIM0, IM_low), IM_high);
        doClip = (IM_low > IM_high);

        scaled = records; n_clipped = 0; s_all = zeros(1,numel(records));
        for kk = 1:numel(records)
            s_raw = targetIM / records(kk).IM;
            if doClip
                s = min(max(s_raw, s_bounds(1)), s_bounds(2));
                if abs(s - s_raw) > 1e-12, n_clipped = n_clipped + 1; end
            else
                s = s_raw;
            end
            s_all(kk) = s;

            scaled(kk).ag    = s * records(kk).ag;
            scaled(kk).PGA   = s * records(kk).PGA;
            scaled(kk).PGV   = s * records(kk).PGV;
            scaled(kk).scale = s;
            scaled(kk).s_clipped = doClip && (abs(s - s_raw) > 1e-12);
            scaled(kk).trimmed   = false;
            scaled(kk).IM = compute_IM_local(scaled(kk).t, scaled(kk).ag, IM_mode, T1, band_fac, band_N);
        end

        err = abs([scaled.IM] - targetIM) / max(targetIM, eps) * 100;
        modeStr = tern_local(strcmpi(IM_mode,'band'),'band','PSA@T1');
        clipCount = n_clipped * doClip;
        if getfield_default_local(opts,'verbose',true)
            fprintf('Hedef IM = %.3f (%s). Maks hata = %.2f%% | uygun aralık=[%.3f, %.3f] | s_min=%.2f s_max=%.2f | KIRPILAN=%d\n', ...
                targetIM, modeStr, max(err), IM_low, IM_high, min(s_all), max(s_all), clipCount);
        end
        meta = struct('IM_mode', IM_mode, 'band_fac', band_fac, 's_bounds', s_bounds);
        if exist('dropped','var'), meta.TRIM_names = dropped; else, meta.TRIM_names = {}; end
    end
end

function IM = compute_IM_local(t, ag, mode, T1, band_fac, band_N)
    zeta = 0.05;
    if strcmpi(mode,'band')
        Tgrid = linspace(band_fac(1)*T1, band_fac(2)*T1, band_N);
        Sa = zeros(size(Tgrid));
        for ii = 1:numel(Tgrid)
            Sa(ii) = calc_psa_local(t, ag, Tgrid(ii), zeta);
        end
        IM = exp(mean(log(Sa + eps)));
    else
        IM = calc_psa_local(t, ag, T1, zeta);
    end
end

function Sa = calc_psa_local(t, ag, T, zeta)
    w = 2*pi / T;
    agf = griddedInterpolant(t, ag, 'linear', 'nearest');
    odef = @(tt, y)[ y(2);
                     -2*zeta*w*y(2) - w*w*y(1) - agf(tt) ];
    y0 = [0;0];
    [~, y] = ode45(odef, t, y0);
    x  = y(:,1);
    xd = y(:,2);
    xdd = -2*zeta*w*xd - w*w*x - ag;
    abs_acc = xdd + ag;
    Sa = max(abs(abs_acc));
end

function s = tern_local(c,a,b)
    if c, s=a; else, s=b; end
end

function [x,a_rel] = solve_linear_mck_local(t, ag, M, C, K)
    nloc = size(M,1);
    r = ones(nloc,1);
    ag_interp = griddedInterpolant(t, ag, 'linear', 'nearest');
    odefun = @(tt,z) [ z(nloc+1:end); ...
        M \ (-C*z(nloc+1:end) - K*z(1:nloc) - M*r*ag_interp(tt)) ...
    ];
    z0 = zeros(2*nloc,1);
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    sol = ode15s(odefun, [t(1) t(end)], z0, opts);
    z = deval(sol, t).';
    x = z(:,1:nloc);
    a_rel = (-(M \ (C*z(:,nloc+1:end).' + K*z(:,1:nloc).')).' - ag.*r.');
end

function [x,a_rel,ts] = mck_with_damper_local(t,ag,M,C,K, k_sd,c_lam0,Lori, orf,rho,Ap,Ao, k_hyd_single,k_p_single, mu_ref, ...
    thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_max,Lgap, cp_oil,cp_steel, ...
    steel_to_oil_mass_ratio, story_mask, n_dampers_per_story, V_oil_per, cfg, dP_max)

    if nargin < 18 || isempty(c_lam_min), c_lam_min = 0; end
    if nargin < 19 || isempty(c_lam_max), c_lam_max = inf; end
    if nargin < 27 || isempty(V_oil_per), V_oil_per = Ap*(2*Lgap); end
    if nargin < 28 || ~isstruct(cfg), cfg = struct(); end
    if nargin < 29 || isempty(dP_max), dP_max = NaN; end

    nloc = size(M,1); r = ones(nloc,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');

    nStories = nloc-1;
    mask = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi = mask .* ndps;
    multi = multi(:);
    Nvec = (1:nStories).'; Mvec = (2:nloc).';

    if isscalar(k_sd)
        k_sd_vec = k_sd * ones(1,nStories);
    else
        k_sd_vec = k_sd(:).';
        if numel(k_sd_vec)~=nStories
            k_sd_vec = k_sd_vec(1) * ones(1,nStories);
        end
    end
    k_sd_vec(~multi.') = 0;

    if numel(V_oil_per)==1
        V_oil_vec = V_oil_per * ones(nStories,1);
    else
        V_oil_vec = V_oil_per(:);
        if numel(V_oil_vec)~=nStories
            V_oil_vec = V_oil_vec(1) * ones(nStories,1);
        end
    end
    if isscalar(k_hyd_single)
        k_hyd_single_vec = k_hyd_single * ones(1,nStories);
    else
        k_hyd_single_vec = k_hyd_single(:).';
        if numel(k_hyd_single_vec)~=nStories
            k_hyd_single_vec = k_hyd_single_vec(1) * ones(1,nStories);
        end
    end
    if isscalar(k_p_single)
        k_p_single_vec = k_p_single * ones(1,nStories);
    else
        k_p_single_vec = k_p_single(:).';
        if numel(k_p_single_vec)~=nStories
            k_p_single_vec = k_p_single_vec(1) * ones(1,nStories);
        end
    end
    active_mask = multi.' > 0;
    k_hyd_single_vec(~active_mask) = 0;
    k_p_single_vec(~active_mask) = 0;

    idx_To_start = 2*nloc + 1;
    idx_Ts_start = 2*nloc + nStories + 1;
    z0 = zeros(2*nloc + 2*nStories,1);
    z0(idx_To_start:idx_To_start+nStories-1) = T0_C;
    z0(idx_Ts_start:idx_Ts_start+nStories-1) = T0_C;

    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);

    c_lam0 = min(max(c_lam0, c_lam_min), c_lam_max);
    c_lam_factor = (128/pi) * Lori * (Ap.^2) / (orf.d_o.^4);

    mu_min_ratio = NaN; mu_max_ratio = NaN;
    if isfield(cfg,'num')
        mu_min_ratio = getfield_default_local(cfg.num,'mu_min_phys',NaN);
        mu_max_ratio = getfield_default_local(cfg.num,'mu_max_phys',NaN);
    end
    mu_floor = mu_ref;
    if isfinite(mu_min_ratio)
        mu_floor = max(mu_ref * mu_min_ratio, 1e-6);
    end
    mu_cap = Inf;
    if isfinite(mu_max_ratio)
        mu_cap = max(mu_ref * mu_max_ratio, mu_floor);
    end

    hA_os    = getfield_default_local(thermal,'hA_os',    thermal.hA_W_perK);
    hA_o_env = getfield_default_local(thermal,'hA_o_env', thermal.hA_W_perK);
    hA_s_env = getfield_default_local(thermal,'hA_s_env', thermal.hA_W_perK);
    T_env    = getfield_default_local(thermal,'T_env_C', T0_C);
    T_floor  = getfield_default_local(thermal,'T_floor_C', -Inf);
    T_ceil   = getfield_default_local(thermal,'T_ceil_C',  Inf);

    function dzdt = damper_ode(tt,z)
        x_loc = z(1:nloc);
        v_loc = z(nloc+1:2*nloc);
        T_o_loc = z(idx_To_start:idx_To_start+nStories-1);
        T_s_loc = z(idx_Ts_start:idx_Ts_start+nStories-1);

        drift_loc = x_loc(Mvec) - x_loc(Nvec);
        dvel_loc  = v_loc(Mvec) - v_loc(Nvec);

        [F_story_loc, F_orf_tot, F_lam_tot, dP_tot_loc, Cd_loc, mu_loc, c_lam_loc, Q_per_loc, P_visc_per_loc, P_orf_per_loc] = ...
            damper_response(drift_loc, dvel_loc, T_o_loc, T_s_loc, tt);

        Fd = zeros(nloc,1);
        Fd(Nvec) = Fd(Nvec) - F_story_loc;
        Fd(Mvec) = Fd(Mvec) + F_story_loc;

        acc = M \ (-C*v_loc - K*x_loc - Fd - M*r*agf(tt));

        [dT_o_dt, dT_s_dt] = temperature_rates(T_o_loc, T_s_loc, P_visc_per_loc, P_orf_per_loc);

        dzdt = [v_loc; acc; dT_o_dt; dT_s_dt];
    end

    sol  = ode15s(@damper_ode,[t(1) t(end)],z0,opts);
    z    = deval(sol,t).';
    x    = z(:,1:nloc);
    v    = z(:,nloc+1:2*nloc);
    T_o  = z(:,idx_To_start:idx_To_start+nStories-1);
    T_s  = z(:,idx_Ts_start:idx_Ts_start+nStories-1);

    drift = x(:,Mvec) - x(:,Nvec);
    dvel  = v(:,Mvec) - v(:,Nvec);

    [F_story, F_orf_tot, F_lam_tot, dP_tot, Cd_hist, mu_hist, c_lam_hist, Q_per, P_visc_per, P_orf_per] = ...
        damper_response(drift, dvel, T_o, T_s, t);

    F = zeros(numel(t),nloc);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a_rel = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    P_total = (P_visc_per + P_orf_per) .* multi.';
    P_orf_tot = P_orf_per .* multi.';
    P_struct = sum(F_story .* dvel, 2);
    E_orf = cumtrapz(t, sum(P_orf_tot,2));
    E_struct = cumtrapz(t, P_struct);

    F_lin_hist = drift .* k_sd_vec;
    F_hyd_hist = drift .* (k_hyd_single_vec .* multi.');
    F_shear_hist = drift .* (k_p_single_vec .* multi.');
    ts = struct('dvel', dvel, 'story_force', F_story, ...
        'F_components', struct('F_orf', F_orf_tot, 'F_lam', F_lam_tot, 'F_lin', F_lin_hist, ...
            'F_hyd', F_hyd_hist, 'F_shear', F_shear_hist), ...
        'Q', Q_per, 'dP_total', dP_tot, 'Cd', Cd_hist, ...
        'mu_inst', mu_hist, 'c_lam_inst', c_lam_hist, ...
        'T_o', T_o, 'T_s', T_s, ...
        'P_visc', P_visc_per, 'P_orf', P_orf_per, ...
        'E_orf', E_orf, 'E_struct', E_struct);

    function [F_story_out, F_orf_tot_out, F_lam_tot_out, dP_tot_out, Cd_out, mu_out, c_lam_out, Q_per_out, P_visc_per_out, P_orf_per_out] = ...
            damper_response(drift_loc, dvel_loc, T_o_loc, T_s_loc, time_loc)

        sz = size(drift_loc);
        drift_vec = reshape(drift_loc, [], nStories);
        dvel_vec  = reshape(dvel_loc,  [], nStories);
        To_vec    = reshape(T_o_loc,   [], nStories);
        Ts_vec    = reshape(T_s_loc,   [], nStories); %#ok<NASGU>

        mu_temp = mu_ref * exp(b_mu * (To_vec - T_ref_C));
        mu_temp = max(mu_temp, mu_floor);
        if isfinite(mu_cap)
            mu_temp = min(mu_temp, mu_cap);
        end
        mu_out = mu_temp;

        c_lam_single_inst = c_lam_factor * mu_temp;
        c_lam_total_inst = c_lam_single_inst .* multi.';
        damper_mask = multi.' > 0;
        if any(damper_mask)
            c_lam_total_inst(:,damper_mask) = min(max(c_lam_total_inst(:,damper_mask), c_lam_min), c_lam_max);
        end
        c_lam_total_inst(:,~damper_mask) = 0;
        c_lam_out = c_lam_total_inst;

        c_lam_single_inst = c_lam_total_inst ./ max(multi.', 1e-12);
        c_lam_single_inst(multi.'==0) = 0;
        F_lam_single = c_lam_single_inst .* dvel_vec;
        F_lam_tot = F_lam_single .* multi.';

        F_hyd_single = drift_vec .* k_hyd_single_vec;
        F_shear_single = drift_vec .* k_p_single_vec;
        F_hyd_tot = F_hyd_single .* multi.';
        F_shear_tot = F_shear_single .* multi.';
        F_lin_tot = F_hyd_tot + F_shear_tot;
        F_lin_per = zeros(size(F_lin_tot));
        mask_multi = multi.' > 0;
        F_lin_per(:,mask_multi) = F_lin_tot(:,mask_multi) ./ multi(mask_multi).';

        params_orf = struct('Ap',Ap,'Ao',Ao,'orf',orf,'rho',rho, ...
            'F_hyd_single',F_hyd_single,'F_visc_single',F_lam_single, ...
            'drift',drift_vec, ...
            'cfgPF',getfield_default_local(cfg,'PF',struct()), ...
            'k_hyd_single',k_hyd_single_vec,'multi',multi.');
        [F_orf_single, dP_tot, Q_per_val, P_orf_single_val, Cd_val] = ...
            calc_orifice_force_local(dvel_vec, params_orf, mu_temp, dP_max);
        F_orf_tot = F_orf_single .* multi.';

        F_resist = F_orf_tot + F_lam_tot;
        F_story_val = F_lin_tot + F_resist;

        P_visc_single = c_lam_single_inst .* (dvel_vec.^2);
        P_orf_single = P_orf_single_val;

        F_story_out = reshape(F_story_val, sz);
        F_orf_tot_out = reshape(F_orf_tot, sz);
        F_lam_tot_out = reshape(F_lam_tot, sz);
        dP_tot_out = reshape(dP_tot, sz);
        Cd_out = reshape(Cd_val, sz);
        mu_out = reshape(mu_temp, sz);
        c_lam_out = reshape(c_lam_total_inst, sz);
        Q_per_out = reshape(Q_per_val, sz);
        P_visc_per_out = reshape(P_visc_single, sz);
        P_orf_per_out = reshape(P_orf_single, sz);
    end

    function [dTo_dt, dTs_dt] = temperature_rates(T_o_loc, T_s_loc, P_visc_per_loc, P_orf_per_loc)
        % Newton tipi ısı transferi: yağ/çelik kütleleri ortam sıcaklığına doğru
        % hA katsayılarıyla soğur; sadece T_floor/T_ceil sayısal güvenlik için
        % kullanılır, fiziksel olarak ortam sıcaklığının altına soğuma mümkündür.
        T_o_loc = T_o_loc(:); T_s_loc = T_s_loc(:);
        P_visc_per_loc = P_visc_per_loc(:);
        P_orf_per_loc = P_orf_per_loc(:);
        dTo_dt = zeros(nStories,1);
        dTs_dt = zeros(nStories,1);
        V_story = V_oil_vec .* multi;
        for jj = 1:nStories
            if multi(jj) <= 0 || V_story(jj) <= 0
                continue;
            end
            C_oil_story = max(rho * V_story(jj) * cp_oil, eps);
            C_steel_story = max(steel_to_oil_mass_ratio * rho * V_story(jj) * cp_steel, eps);
            power_in = (P_visc_per_loc(jj) + P_orf_per_loc(jj)) * multi(jj);
            dT_o = ( power_in - hA_os*(T_o_loc(jj) - T_s_loc(jj)) - hA_o_env*(T_o_loc(jj) - T_env) ) / C_oil_story;
            dT_s = ( +hA_os*(T_o_loc(jj) - T_s_loc(jj)) - hA_s_env*(T_s_loc(jj) - T_env) ) / C_steel_story;
            if T_o_loc(jj) <= T_floor && dT_o < 0, dT_o = 0; end
            if T_o_loc(jj) >= T_ceil  && dT_o > 0, dT_o = 0; end
            if T_s_loc(jj) <= T_floor && dT_s < 0, dT_s = 0; end
            if T_s_loc(jj) >= T_ceil  && dT_s > 0, dT_s = 0; end
            dTo_dt(jj) = dT_o;
            dTs_dt(jj) = dT_s;
        end
    end
end

function [F_orf, dP_total, Q, P_orf_per, Cd] = calc_orifice_force_local(dvel, params, mu, dP_max)
    % Fiziksel hazne modeli: piston akışı Q, lineer hidrolik yay, viskoz sürtünme
    % ve gazlı akümülatör (cfg.PF) ile tanımlanan arka basınç kombinasyonu.
    % Gaz basıncı, lineer rezervuar katsayıları ve bleed yolundaki direnç
    % paralel bağlanmış elemanlar olarak çözülür; böylece gaz precharge ve
    % rezervuar sertlikleri arka basıncı (p_back) ambient üzerine taşıyabilir.
    if nargin < 3 || isempty(mu)
        mu = params.mu;
    end
    if nargin < 4
        dP_max = NaN;
    end
    multi_raw = getfield_default_local(params,'multi',ones(1,size(dvel,2)));
    active_mask_row = multi_raw > 0;
    if size(active_mask_row,1) > 1
        active_mask_row = active_mask_row(1,:);
    end
    multi_mask = repmat(active_mask_row, size(dvel,1), 1);
    Q = params.Ap .* dvel .* multi_mask;
    Ao_eff = max(params.Ao, 1e-12);
    Re = (params.rho .* abs(Q) .* params.orf.d_o) ./ max(Ao_eff .* mu, 1e-12);
    Cd0   = params.orf.Cd0;
    CdInf = params.orf.CdInf;
    p_exp = params.orf.p_exp;
    Rec   = params.orf.Rec;
    Cd    = CdInf - (CdInf - Cd0) ./ (1 + (Re./max(Rec,1)).^p_exp);
    Cd    = max(min(Cd, 1.2), 0.2);

    dP_kv = 0.5 * params.rho .* (Q ./ max(Cd .* Ao_eff, 1e-12)).^2;

    F_hyd = getfield_default_local(params,'F_hyd_single',zeros(size(Q)));
    F_visc = getfield_default_local(params,'F_visc_single',zeros(size(Q)));
    p_up = params.orf.p_amb + (F_hyd + F_visc) ./ max(params.Ap, 1e-12);

    drift = getfield_default_local(params,'drift',zeros(size(Q)));
    cfgPF = getfield_default_local(params,'cfgPF',struct());
    resv = getfield_default_local(params.orf,'reservoir',struct());
    V_ref = getfield_default_local(resv,'V_ref', getfield_default_local(cfgPF,'V_gas_ref',1e-3));
    p_min = getfield_default_local(resv,'p_min', params.orf.p_cav_eff);
    k_lin = getfield_default_local(resv,'k_lin', 0);
    c_lin = getfield_default_local(resv,'c_lin', 0);
    gamma = getfield_default_local(cfgPF,'gas_gamma',1.2);
    V_gas0 = getfield_default_local(cfgPF,'V_gas_ref', max(V_ref,1e-6));
    p_pre = getfield_default_local(cfgPF,'precharge_Pa', params.orf.p_amb);
    bleed_coeff = getfield_default_local(cfgPF,'bleed_coeff',0);

    deltaV = drift .* params.Ap .* multi_mask;
    V_gas = max(V_gas0 + deltaV, 1e-6);
    % Gaz akümülatörü precharge (cfg.PF.precharge_Pa) artarken p_gas büyür ve
    % hacimsel sertlik kappa_gas üzerinden p_back'i yukarı çeker. Rezervuarın
    % lineer yay/damper katsayıları (k_lin, c_lin) aynı şekilde p_res_lin'i ve
    % onun ağırlığını belirler. Bleed hattı ise bleed_coeff azaldıkça ambiente
    % daha güçlü bağlanarak p_back'i p_min'e doğru çeker.
    p_gas = p_pre .* (max(V_gas0,1e-6) ./ V_gas) .^ gamma;
    p_res_lin = p_min + k_lin .* (deltaV ./ max(V_ref, 1e-6)) + c_lin .* (Q ./ max(params.Ap, 1e-12));
    p_bleed = p_min + bleed_coeff .* (Q ./ max(params.Ao, 1e-12));

    kappa_gas = gamma .* p_gas ./ max(V_gas, 1e-9);
    kappa_lin = max(k_lin ./ max(V_ref, 1e-9), 0) + max(abs(c_lin) ./ max(params.Ap, 1e-12), 0);
    bleed_scale = max(abs(bleed_coeff), 1e-12);
    w_bleed = max(params.Ao, 1e-12) ./ bleed_scale;
    w_bleed(bleed_coeff == 0) = 0;
    w_base = 1e-9;
    w_sum = kappa_gas + kappa_lin + w_bleed + w_base;
    p_back = (kappa_gas .* p_gas + kappa_lin .* p_res_lin + w_bleed .* p_bleed + w_base .* p_min) ./ w_sum;
    p_back = max(p_back, p_min);

    dP_supply = max(p_up - p_back, 0);
    p_cav_floor = max(params.orf.p_cav_eff, 0);
    p_supply_cap = getfield_default_local(params.orf,'p_supply_cap', params.orf.p_amb);
    p_supply_cap = max(p_supply_cap, p_cav_floor);
    p_supply_cap = p_supply_cap + zeros(size(p_up));
    p_up_cav = min(max(p_up, p_cav_floor), p_supply_cap);
    dP_cav = max((p_up_cav - p_cav_floor) .* params.orf.cav_sf, 0);
    dP_allow = max(min(dP_supply, dP_cav), 0);
    if isfinite(dP_max)
        dP_kv = min(dP_kv, dP_max);
        dP_allow = min(dP_allow, dP_max);
    end
    eps_spring = abs(k_lin .* (deltaV ./ max(V_ref, 1e-6)));
    eps_damp = abs(c_lin .* (Q ./ max(params.Ap, 1e-12)));
    eps_bleed = abs(bleed_coeff .* (Q ./ max(params.Ao, 1e-12)));
    eps_base = getfield_default_local(params.orf,'smooth_eps0',1e3);
    epsm = max(eps_base, 0.05 * dP_allow + eps_spring + eps_damp + eps_bleed);
    dP_total = softmin_local(dP_kv, dP_allow, epsm);
    sgn = sign(Q);
    F_orf = dP_total .* params.Ap .* sgn;
    P_orf_per = dP_total .* Q;
end
