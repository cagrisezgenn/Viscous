function x = trial_x_all_from_structs(geom, sh, orf, hyd, therm, num, ga)
%TRIAL_X_ALL_FROM_STRUCTS Build design vector snapshot from structs.
% Layout:
%  set1 (1..6):  [d_o, n_orf, Lori, Dp, Lgap, mu_ref]
%  set2 (7..11): [Cd0, CdInf, Rec, p_exp, PF_gain]
%  set3 (12..15):[Kd, beta0, K_leak, Vmin_fac]
%  set4 (16..19):[d_w, D_m, n_turn, n_parallel]
%
% Inputs:
%   geom, sh, orf, hyd, therm - parameter structures
%   num  - kept for compatibility (unused)
%   ga   - GA configuration structure (for PF_gain)
%
% Output:
%   x  - 19-element design vector assembled from the structs.
%
    if nargin < 7
        ga = struct();
    end
    %#ok<INUSD> % num is unused but kept for compatibility

    x = nan(1,19);
    gfd = @getfield_default;

    % --- set 1 ---
    x(1) = gfd(geom, 'd_o',   2.2e-3);
    x(2) = gfd(orf,  'n_orf', 2);
    x(3) = gfd(geom, 'Lori',  0.03);
    x(4) = gfd(geom, 'Dp',    0.12);
    x(5) = gfd(geom, 'Lgap',  0.20);
    x(6) = gfd(therm,'mu_ref',1.2);

    % --- set 2 ---
    x(7)  = gfd(orf, 'Cd0',    0.6);
    x(8)  = gfd(orf, 'CdInf',  0.9);
    x(9)  = gfd(orf, 'Rec',    3800);
    x(10) = gfd(orf, 'p_exp',  1.1);
    x(11) = gfd(ga,  'PF_gain',1.0);

    % --- set 3 ---
    x(12) = gfd(geom,'Kd',       1.8e9);
    x(13) = gfd(therm,'beta0',   1.6e9);
    x(14) = gfd(hyd, 'K_leak',   0);
    x(15) = gfd(hyd, 'Vmin_fac', 0.97);

    % --- set 4 ---
    x(16) = gfd(sh,  'd_w',       0.018);
    x(17) = gfd(sh,  'D_m',       0.08);
    x(18) = gfd(sh,  'n_turn',    26);
    x(19) = gfd(hyd, 'n_parallel',1);
end
