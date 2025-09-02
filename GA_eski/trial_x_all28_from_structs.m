function x28 = trial_x_all28_from_structs(geom, sh, orf, hyd, therm, num)
%TRIAL_X_ALL28_FROM_STRUCTS Build 28-length design vector snapshot from structs.
% Layout:
%  set1 (1..9):  [d_o, Lori, mu_ref, Kd, d_w, D_m, n_turn, Dp, Lgap]
%  set2 (10..18):[n_orf, Cd0, CdInf, Rec, p_exp, cav_sf, Lh, K_leak, resFactor]
%  set3 (19..28):[n_orf, Cd0, mu_ref, b_mu, beta0, b_beta, hA_os, dP_cap, Vmin_fac, resFactor]

    x28 = nan(1,28);
    gfd = @getfield_default;
    % --- set 1 ---
    x28(1) = gfd(geom,'d_o',   2.2e-3);
    x28(2) = gfd(geom,'Lori',  0.03);
    x28(3) = gfd(therm,'mu_ref', 1.2);
    x28(4) = gfd(geom,'Kd',    1.8e9);
    x28(5) = gfd(sh,  'd_w',   0.018);
    x28(6) = gfd(sh,  'D_m',   0.08);
    x28(7) = gfd(sh,  'n_turn', 26);
    x28(8) = gfd(geom,'Dp',    0.12);
    x28(9) = gfd(geom,'Lgap',  0.20);

    % --- set 2 ---
    x28(10) = gfd(orf, 'n_orf', 2);
    x28(11) = gfd(orf, 'Cd0',   0.6);
    x28(12) = gfd(orf, 'CdInf', 0.9);
    x28(13) = gfd(orf, 'Rec',   3800);
    x28(14) = gfd(orf, 'p_exp', 1.1);
    x28(15) = gfd(orf, 'cav_sf',1.20);
    x28(16) = gfd(hyd, 'Lh',    9e-4);
    x28(17) = gfd(hyd, 'K_leak',0);
    x28(18) = gfd(therm,'resFactor',22);

    % --- set 3 ---
    x28(19) = gfd(orf,  'n_orf', 2);
    x28(20) = gfd(orf,  'Cd0',   0.6);
    x28(21) = gfd(therm,'mu_ref',1.2);
    x28(22) = gfd(therm,'b_mu',  -0.011);
    x28(23) = gfd(therm,'beta0', 1.6e9);
    x28(24) = gfd(therm,'b_beta',-0.0035);
    x28(25) = gfd(therm,'hA_os', 800);
    x28(26) = gfd(num,  'dP_cap',5e7);
    x28(27) = gfd(hyd,  'Vmin_fac',0.97);
    x28(28) = gfd(therm,'resFactor',22);
end

