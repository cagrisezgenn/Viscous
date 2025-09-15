%% SDOF analysis using GA best parameters
% Runs a single simulation without optimization using the best design
% parameters (xbest) stored in ga_front.csv. Produces displacement and
% acceleration histories for the 10th floor along with peak metrics.

clear; clc; close all;

%% Load base parameters
parametreler;  % loads structural and damper parameters and computes T1

%% Override with GA best parameters from GA output
try
    if exist('ga_knee.csv','file')
        tbl = readtable('ga_knee.csv');
    else
        tbl = readtable('ga_front.csv');
    end
    xb = tbl(1,:);

    d_o   = xb.d_o_mm/1000;      % [m]
    n_orf = xb.n_orf;
    cfg.PF.tau  = xb.PF_tau;
    cfg.PF.gain = xb.PF_gain;
    cfg.PF.t_on = xb.PF_t_on;
    orf.Cd0  = xb.Cd0;
    orf.CdInf= xb.CdInf;
    orf.p_exp= xb.p_exp;
    Lori  = xb.Lori_mm/1000;     % [m]
    thermal.hA_W_perK = xb.hA_W_perK;
    Dp    = xb.Dp_mm/1000;       % [m]
    d_w   = xb.d_w_mm/1000;      % [m]
    D_m   = xb.D_m_mm/1000;      % [m]
    n_turn= xb.n_turn;
    mu_ref= xb.mu_ref;

    % Recompute derived damper constants
    params_tmp = struct('Dp',Dp,'d_w',d_w,'D_m',D_m,'n_turn',n_turn, ...
        'mu_ref',mu_ref,'Lori',Lori,'Lgap',Lgap,'Kd',Kd,'Ebody',Ebody, ...
        'Gsh',Gsh,'orf',struct('d_o',d_o));
    params_tmp = build_params(params_tmp);
    Ap    = params_tmp.Ap;
    k_sd  = params_tmp.k_sd;
    c_lam0= params_tmp.c_lam0;
    orf.d_o = params_tmp.orf.d_o;
    A_o   = n_orf * (pi*orf.d_o^2/4);
    Qcap_big = max(orf.CdInf*A_o,1e-9) * sqrt(2*1.0e9/rho);
    c_lam_min = max(c_lam_min_abs, c_lam_min_frac*c_lam0);
catch ME
    warning('Unable to apply GA parameters: %s', ME.message);
end

%% Ground motion input (scaled by T1)
[recs_raw, recs] = load_ground_motions(T1);
rec = recs(1);  t = rec.t; ag = rec.ag;
win = Utils.make_arias_window(t, ag); t5 = win.t5; t95 = win.t95;

%% Damperless and damper responses
[x0,a_rel0] = Utils.lin_MCK(t, ag, M, C0, K);
[x_d,a_d,diag_d] = mck_with_damper(t, ag, M, C0, K, k_sd, c_lam0, Lori, ...
    orf, rho, Ap, A_o, Qcap_big, mu_ref, thermal, T0_C, T_ref_C, b_mu, ...
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
