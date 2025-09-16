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
story_mask = params.story_mask;  resFactor = params.resFactor;
Kd = params.Kd; Ebody = params.Ebody; Gsh = params.Gsh;

Ap = params.Ap;   Ao = params.Ao;   k_sd = params.k_sd;
c_lam0 = params.c_lam0;  c_lam_cap = params.c_lam_cap;
c_lam_min_abs = params.c_lam_min_abs;  c_lam_min_frac = params.c_lam_min_frac;
c_lam_min = params.c_lam_min;
Qcap_big = params.Qcap_big;

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
        'c_lam_min_abs',c_lam_min_abs,'c_lam_min_frac',c_lam_min_frac);

    params_override = build_params_local(params_override);

    Ap = params_override.Ap;
    Ao = params_override.Ao;
    k_sd = params_override.k_sd;
    c_lam0 = params_override.c_lam0;
    c_lam_min = params_override.c_lam_min;
    Qcap_big = params_override.Qcap_big;
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
    orf, rho, Ap, Ao, Qcap_big, mu_ref, thermal, T0_C, T_ref_C, b_mu, ...
    c_lam_min, c_lam_cap, Lgap, cp_oil, cp_steel, steel_to_oil_mass_ratio, ...
    story_mask, n_dampers_per_story, resFactor, cfg);

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
for i = 1:nStories
    idx = [i,i+1];
    k_eq = k_sd * multi(i);
    c_eq = diag_d.c_lam * multi(i);
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
    params.orf = struct('Cd0',0.61,'CdInf',0.80,'Rec',3000,'p_exp',1.1, ...
        'p_amb',1.0e5,'p_cav_eff',2.0e3,'cav_sf',0.90,'d_o',3.0e-3,'veps',0.10);
    params.Lori = 0.10;
    params.mu_ref = 0.9;

    params.Kd = 1.6e9;
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
        'tol_K',0.5,'relax',0.5,'dT_max',80);

    params.steel_to_oil_mass_ratio = 1.5;
    params.n_dampers_per_story = 1;
    params.story_mask = ones(n-1,1);
    params.cp_oil = 1800;
    params.cp_steel = 500;
    params.resFactor = 12;

    params.c_lam_cap = 2e7;
    params.c_lam_min_frac = 0.05;
    params.c_lam_min_abs = 1e5;

    cfg = struct();
    cfg.PF = struct('mode','ramp','tau',1.0,'gain',0.85,'t_on',0,'auto_t_on',true,'k',0.01);
    cfg.on = struct('pressure_force',true,'mu_floor',false);
    cfg.compat_simple = false;
    cfg.num = struct('softmin_eps',1e4,'mu_min_phys',0.6,'dP_cap',NaN);
    params.cfg = cfg;

    params = build_params_local(params);
end

function params = build_params_local(params)
    params = recompute_damper_params_local(params);
    if isfield(params,'orf') && isfield(params.orf,'CdInf') && ...
            isfield(params,'Ao') && isfield(params,'rho')
        params.Qcap_big = max(params.orf.CdInf * params.Ao, 1e-9) * ...
            sqrt(2*1.0e9 / params.rho);
    else
        params.Qcap_big = 0;
    end
    if isfield(params,'c_lam_min_abs') && isfield(params,'c_lam_min_frac') && ...
            isfield(params,'c_lam0')
        params.c_lam_min = max(params.c_lam_min_abs, ...
            params.c_lam_min_frac * params.c_lam0);
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

    Ap = pi * params.Dp^2 / 4;
    Ao_single = pi * params.orf.d_o^2 / 4;
    Ao = params.n_orf * Ao_single;
    Ap_eff = nd * Ap;
    Ao_eff = nd * Ao;

    rho_loc = getfield_default_local(params,'rho',850);
    Lh = rho_loc * params.Lori / max(Ao^2, 1e-18);

    k_h = params.Kd * Ap^2 / params.Lgap;
    k_s = params.Ebody * Ap / params.Lgap;
    k_hyd = 1 / (1/k_h + 1/k_s);
    k_p = params.Gsh * params.d_w^4 / (8 * params.n_turn * params.D_m^3);
    k_sd_simple = k_hyd + k_p;
    k_sd_adv    = nd * (k_hyd + k_p);

    c_lam0 = 12 * params.mu_ref * params.Lori * Ap^2 / (params.orf.d_o^4);

    params.Ap = Ap;
    params.Ao = Ao;
    params.Ap_eff = Ap_eff;
    params.Ao_eff = Ao_eff;
    params.Lh = Lh;
    params.k_p = k_p;
    params.k_sd_simple = k_sd_simple;
    params.k_sd_adv = k_sd_adv;
    params.k_sd = k_sd_adv;
    params.c_lam0 = c_lam0;
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

function w = pf_weight_local(t, cfg)
    if nargin < 2 || ~isstruct(cfg), cfg = struct(); end
    if ~isfield(cfg,'on') || ~isstruct(cfg.on), cfg.on = struct(); end
    if ~isfield(cfg.on,'pressure_force'), cfg.on.pressure_force = true; end
    if ~isfield(cfg,'PF') || ~isstruct(cfg.PF), cfg.PF = struct(); end
    if ~isfield(cfg.PF,'t_on'), cfg.PF.t_on = 0; end
    if ~isfield(cfg.PF,'tau'),  cfg.PF.tau  = 1.0; end
    if ~isfield(cfg,'compat_simple'), cfg.compat_simple = true; end
    k = getfield_default_local(cfg.PF,'k',0.01);
    tau_floor = 1e-6;
    t = double(t);
    if cfg.compat_simple
        dt  = max(t - cfg.PF.t_on, 0);
        tau = max(cfg.PF.tau, tau_floor);
        w_local = 1 - exp(-dt ./ tau);
    else
        k = max(k, tau_floor);
        sp_dt = (log1p(exp(-abs((t - cfg.PF.t_on)./k))) + max((t - cfg.PF.t_on)./k, 0));
        dt  = sp_dt .* k;
        sp_tau = (log1p(exp(-abs((cfg.PF.tau - tau_floor)./k))) + max((cfg.PF.tau - tau_floor)./k, 0));
        tau = sp_tau .* k + tau_floor;
        w_local = 1 - exp(-dt ./ tau);
    end
    w = cfg.on.pressure_force .* w_local;
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

function [x,a_rel,ts] = mck_with_damper_local(t,ag,M,C,K, k_sd,c_lam0,Lori, orf,rho,Ap,Ao,Qcap, mu_ref, ...
    thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, story_mask, ...
    n_dampers_per_story, resFactor, cfg)

    nloc = size(M,1); r = ones(nloc,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0 = zeros(2*nloc,1);
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    nStories = nloc-1;
    mask = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi= (mask .* ndps).';
    Nvec = 1:nStories; Mvec = 2:nloc;

    Tser = T0_C*ones(numel(t),1);
    mu_abs = mu_ref;
    c_lam = min(max(c_lam0, c_lam_min), c_lam_cap);

    odef = @(tt,z) [ z(nloc+1:end); M \ ( -C*z(nloc+1:end) - K*z(1:nloc) - dev_force(tt,z(1:nloc),z(nloc+1:end),c_lam,mu_abs) - M*r*agf(tt) ) ];
    sol  = ode15s(odef,[t(1) t(end)],z0,opts);
    z    = deval(sol,t).';
    x    = z(:,1:nloc); v = z(:,nloc+1:end);

    drift = x(:,Mvec) - x(:,Nvec);
    dvel  = v(:,Mvec) - v(:,Nvec);
    F_lin = k_sd*drift;

    Qcap_eff = Qcap;
    if isfield(cfg,'num') && isfield(cfg.num,'Qcap_scale') && isfinite(cfg.num.Qcap_scale)
        Qcap_eff = max(1e-9, Qcap * cfg.num.Qcap_scale);
    end
    orf_loc = orf;
    if isfield(cfg,'num') && isfield(cfg.num,'softmin_eps') && isfinite(cfg.num.softmin_eps)
        orf_loc.softmin_eps = cfg.num.softmin_eps;
    end
    params_orf = struct('Ap',Ap,'Qcap',Qcap_eff,'orf',orf_loc,'rho',rho,...
                        'Ao',Ao,'mu',mu_abs,'F_lin',F_lin,'Lori',Lori);
    [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force_local(dvel, params_orf);

    qmag_loc = Qcap_eff * tanh( (Ap/Qcap_eff) * sqrt(dvel.^2 + orf.veps^2) );
    Re_loc   = (rho .* qmag_loc .* max(orf.d_o,1e-9)) ./ max(Ao*mu_abs,1e-9);
    Cd_loc0  = orf.Cd0; Cd_locInf = orf.CdInf; Rec_loc = orf.Rec; pexp_loc = orf.p_exp;
    Cd_loc = Cd_locInf - (Cd_locInf - Cd_loc0) ./ (1 + (Re_loc./max(Rec_loc,1)).^pexp_loc);
    Cd_loc = max(min(Cd_loc,1.2),0.2);
    dP_kv_loc = 0.5*rho .* ( qmag_loc ./ max(Cd_loc.*Ao,1e-12) ).^2;
    p_up_loc  = orf.p_amb + abs(F_lin)./max(Ap,1e-12);
    dP_cav_loc= max( (p_up_loc - orf.p_cav_eff).*orf.cav_sf, 0 );
    F_p = F_lin + F_orf;

    dp_pf = (c_lam*dvel + (F_p - k_sd*drift)) ./ Ap;
    if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
        s = tanh(20*dvel);
        dp_pf = s .* max(0, s .* dp_pf);
    end
    w_pf_vec = pf_weight_local(t, cfg) * cfg.PF.gain;
    F_p = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;

    F_story = F_p;
    P_visc_per = c_lam * (dvel.^2);
    P_sum = sum( (P_visc_per + P_orf_per) .* multi, 2 );
    P_orf_tot = sum(P_orf_per .* multi, 2);
    P_struct_tot = sum(F_story .* dvel, 2);
    E_orf = cumtrapz(t, P_orf_tot);
    E_struct = cumtrapz(t, P_struct_tot);

    nDtot = sum(multi);
    V_oil_per = resFactor*(Ap*(2*Lgap));
    m_oil_tot = nDtot*(rho*V_oil_per);
    m_steel_tot = steel_to_oil_mass_ratio*m_oil_tot;
    C_oil   = max(m_oil_tot*cp_oil,   eps);
    C_steel = max(m_steel_tot*cp_steel, eps);
    T_o = Tser; T_s = T0_C*ones(numel(t),1);
    hA_os   = getfield_default_local(thermal, 'hA_os',    thermal.hA_W_perK);
    hA_o_env= getfield_default_local(thermal, 'hA_o_env', thermal.hA_W_perK);
    hA_s_env= getfield_default_local(thermal, 'hA_s_env', thermal.hA_W_perK);
    dtv = diff(t);
    for kk=1:numel(t)-1
        Pk = 0.5*(P_sum(kk)+P_sum(kk+1));
        dT_o = ( Pk - hA_os*(T_o(kk)-T_s(kk)) - hA_o_env*(T_o(kk)-thermal.T_env_C) ) / C_oil;
        dT_s = ( + hA_os*(T_o(kk)-T_s(kk)) - hA_s_env*(T_s(kk)-thermal.T_env_C) ) / C_steel;
        T_o(kk+1) = T_o(kk) + dtv(kk)*dT_o;
        T_s(kk+1) = T_s(kk) + dtv(kk)*dT_s;
        T_o(kk+1) = min(max(T_o(kk+1), T0_C), T0_C + thermal.dT_max);
        T_s(kk+1) = min(max(T_s(kk+1), T0_C), T0_C + thermal.dT_max);
    end
    mu = mu_ref*exp(b_mu*(T_o - T_ref_C));

    F = zeros(numel(t),nloc);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a_rel = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    ts = struct('dvel', dvel, 'story_force', F_story, 'Q', Q, ...
        'dP_orf', dP_orf, 'PF', F_p, 'cav_mask', dP_orf < 0, 'P_sum', P_sum, ...
        'E_orf', E_orf, 'E_struct', E_struct, 'T_oil', T_o, 'mu', mu, 'c_lam', c_lam);

    function Fd = dev_force(tt,x_,v_,c_lam_loc,mu_abs_loc)
        drift_ = x_(Mvec) - x_(Nvec);
        dvel_  = v_(Mvec) - v_(Nvec);
        F_lin_ = k_sd*drift_;
        params_loc = struct('Ap',Ap,'Qcap',Qcap,'orf',orf,'rho',rho,...
                        'Ao',Ao,'mu',mu_abs_loc,'F_lin',F_lin_,'Lori',Lori);
        [F_orf_, ~, ~, ~] = calc_orifice_force_local(dvel_, params_loc);
        dp_pf_ = (c_lam_loc*dvel_ + F_orf_) ./ Ap;
        if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
            s = tanh(20*dvel_);
            dp_pf_ = s .* max(0, s .* dp_pf_);
        end
        w_pf = pf_weight_local(tt,cfg) * cfg.PF.gain;
        F_p_ = k_sd*drift_ + (w_pf .* dp_pf_) * Ap;
        F_story_ = F_p_;
        Fd = zeros(nloc,1);
        Fd(Nvec) = Fd(Nvec) - F_story_;
        Fd(Mvec) = Fd(Mvec) + F_story_;
    end
end

function [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force_local(dvel, params)
    qmag = params.Qcap * tanh( (params.Ap/params.Qcap) * sqrt(dvel.^2 + params.orf.veps^2) );
    Re   = (params.rho .* qmag .* max(params.orf.d_o,1e-9)) ./ max(params.Ao*params.mu,1e-9);
    Cd0   = params.orf.Cd0;
    CdInf = params.orf.CdInf;
    p_exp = params.orf.p_exp;
    Rec   = params.orf.Rec;
    Cd    = CdInf - (CdInf - Cd0) ./ (1 + (Re./max(Rec,1)).^p_exp);
    Cd    = max(min(Cd, 1.2), 0.2);

    dP_kv  = 0.5*params.rho .* ( qmag ./ max(Cd.*params.Ao,1e-12) ).^2;
    p_up   = params.orf.p_amb + abs(params.F_lin)./max(params.Ap,1e-12);
    dP_cav = max( (p_up - params.orf.p_cav_eff).*params.orf.cav_sf, 0 );
    epsm = 1e5;
    if isfield(params,'orf') && isfield(params.orf,'softmin_eps') && isfinite(params.orf.softmin_eps)
        epsm = params.orf.softmin_eps;
    end
    dP_orf = softmin_local(dP_kv, dP_cav, epsm);
    sgn = dvel ./ sqrt(dvel.^2 + params.orf.veps^2);
    F_orf = dP_orf .* params.Ap .* sgn;

    Q = qmag;
    P_orf_per = dP_kv .* qmag;
end