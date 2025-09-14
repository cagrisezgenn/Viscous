%% SDOF analysis using GA best parameters
% Synced with latest upstream
% Runs a single simulation without optimization using the best design
% parameters (xbest) stored in ga_front.csv. Produces displacement and
% acceleration histories for the 10th floor along with peak metrics.

clear; clc; close all;

%% Load base parameters
parametreler;  % loads structural and damper parameters and computes T1

%% Override with GA best parameters from ga_front.csv
try
    tbl = readtable('ga_front.csv');
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
    params_tmp = Utils.recompute_damper_params(params_tmp);
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
% Call parent mck_with_damper to avoid local "dP_kv_loc" bug
curdir = pwd; cd('..');
[x_lin,a_lin,diag_lin] = mck_with_damper(t, ag, M, C0, K, k_sd, c_lam0, Lori, false, ...
    orf, rho, Ap, A_o, Qcap_big, mu_ref, false, thermal, T0_C, T_ref_C, b_mu, ...
    c_lam_min, c_lam_cap, Lgap, cp_oil, cp_steel, steel_to_oil_mass_ratio, toggle_gain, ...
    story_mask, n_dampers_per_story, resFactor, cfg);
cd(curdir); cd('..');
[x_orf,a_orf,diag_orf] = mck_with_damper(t, ag, M, C0, K, k_sd, c_lam0, Lori, true, ...
    orf, rho, Ap, A_o, Qcap_big, mu_ref, true, thermal, T0_C, T_ref_C, b_mu, ...
    c_lam_min, c_lam_cap, Lgap, cp_oil, cp_steel, steel_to_oil_mass_ratio, toggle_gain, ...
    story_mask, n_dampers_per_story, resFactor, cfg);
cd(curdir);

x10_0   = x0(:,10);  x10_lin = x_lin(:,10);  x10_orf = x_orf(:,10);
a10_0   = a_rel0(:,10) + ag;
a10_lin = a_lin(:,10) + ag;
a10_orf = a_orf(:,10) + ag;

%% Assemble equivalent damping/stiffness matrices for modal check
nStories = n - 1;
mask = story_mask(:); if numel(mask)==1, mask = mask*ones(nStories,1); end
ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
multi = (mask .* ndps);
Kadd = zeros(n); Cl_add = zeros(n); Co_add = zeros(n);
for i = 1:nStories
    idx = [i,i+1];
    k_eq   = k_sd * multi(i);
    c_eq_l = diag_lin.c_lam * multi(i);
    c_eq_o = diag_orf.c_lam * multi(i);
    kM  = k_eq  * [1 -1; -1 1];
    cMl = c_eq_l* [1 -1; -1 1];
    cMo = c_eq_o* [1 -1; -1 1];
    Kadd(idx,idx)  = Kadd(idx,idx)  + kM;
    Cl_add(idx,idx)= Cl_add(idx,idx)+ cMl;
    Co_add(idx,idx)= Co_add(idx,idx)+ cMo;
end
K_tot = K + Kadd;
C_lin = C0 + Cl_add;
C_orf = C0 + Co_add;
[V,D] = eig(K_tot,M); [w2,ord] = sort(diag(D),'ascend');
phi1 = V(:,ord(1)); w1 = sqrt(w2(1));
normM = phi1.' * M * phi1;

%% Plots: 10th floor displacement and acceleration, plus IDR
figure('Name','10. Kat yer değiştirme — ham ivme (ODE-only)','Color','w');
plot(t, x10_0 ,'k','LineWidth',1.4); hold on;
plot(t, x10_lin,'b','LineWidth',1.1);
plot(t, x10_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('x10(t) [m]');
title(sprintf('10-Kat | T1=%.3f s | Arias [%.3f, %.3f] s', T1, t5, t95));
legend('Dampersiz','Lineer damper','Orifisli damper','Location','best');

figure('Name','10. Kat mutlak ivme','Color','w');
plot(t, a10_0 ,'k','LineWidth',1.4); hold on;
plot(t, a10_lin,'b','LineWidth',1.1);
plot(t, a10_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('a10abs(t) [m/s^2]');
legend('Dampersiz','Lineer damper','Orifisli damper','Location','best');

drift0    = x0(:,2:end)   - x0(:,1:end-1);
drift_lin = x_lin(:,2:end) - x_lin(:,1:end-1);
drift_orf = x_orf(:,2:end) - x_orf(:,1:end-1);
IDR0      = max(abs(drift0))./story_height;
IDR_lin   = max(abs(drift_lin))./story_height;
IDR_orf   = max(abs(drift_orf))./story_height;
story_ids = 1:(n-1);
figure('Name','Maksimum IDR','Color','w');
plot(story_ids, IDR0,'k-o','LineWidth',1.4); hold on;
plot(story_ids, IDR_lin,'b-s','LineWidth',1.1);
plot(story_ids, IDR_orf,'r-d','LineWidth',1.0);
grid on; xlabel('Kat'); ylabel('Maks IDR [Delta x/h]');
legend('Dampersiz','Lineer damper','Orifisli damper','Location','best');

%% Summary printout
zeta0 = (phi1.' * C0 * phi1) / (2*w1*normM);
zeta_d = (phi1.' * (C0 + Co_add) * phi1) / (2*w1*normM);

x10_max_0    = max(abs(x10_0));
x10_max_d    = max(abs(x10_orf));
a10abs_max_0 = max(abs(a10_0));
a10abs_max_d = max(abs(a10_orf));

fprintf('Self-check zeta1: %.3f %% (dampersiz) vs %.3f %% (damperli)\n', 100*zeta0, 100*zeta_d);
fprintf('x10_max  (dampersiz)   = %.4g m\n', x10_max_0);
fprintf('x10_max  (damperli)    = %.4g m\n', x10_max_d);
fprintf('a10abs_max  (dampersiz)= %.4g m/s^2\n', a10abs_max_0);
fprintf('a10abs_max  (damperli) = %.4g m/s^2\n', a10abs_max_d);
