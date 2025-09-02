function [lb,ub,int_idx,names] = ga_get_bounds(set_id)
% Yeni 4Ã—8 GA set tanÄ±mÄ± (toplam 32 deÄŸiÅŸken)
% Set 1 â€” Orifis + AkÄ±ÅŸ: [Cd0, CdInf, Rec, p_exp, cav_sf, K_leak, Vmin_fac, d_o, Lori]
% Set 2 â€” Geometri + Hidrolik: [Dp, Lgap, Kd, dP_cap, F_cap, n_turn, d_w, D_m, mu_ref, beta0, resFactor]
% Set 3 â€” Termal Ã–z: [b_mu, b_beta, rho_ref, alpha_rho, T_ref_C, antoine_A, antoine_B, antoine_C]
% Set 4 â€” Termal SÄ±nÄ±r/IC + IsÄ± GeÃ§iÅŸi: [hA_os, hA_o_env, hA_s_env, cp_oil, cp_steel, T_env_C, T0_C, Ts0_C]

    switch double(set_id)
        case 1
            % [Cd0, CdInf, Rec, p_exp, cav_sf, K_leak, Vmin_fac, d_o, Lori]
            lb = [ 0.30,  0.60, 2533.333333, 0.55, 0.525, 1.0e-9, 0.93, 1.0e-4, 1.0e-4];
            ub = [ 0.90,  1.35, 5700.000000, 1.65, 1.575, 6.0e-9, 0.98, 5.0e-3, 1.0e-2];
            int_idx = [];
            names = {'Cd0','CdInf','Rec','p_exp','cav_sf','K_leak','Vmin_fac','d_o','Lori'};

        case 2
            % [Dp, Lgap, Kd, dP_cap, F_cap, n_turn, d_w, D_m, mu_ref, beta0, resFactor]
            lb = [0.076, 0.0615, 1.08e9, 3.333333333e7, 1e6, 15, 0.01, 0.05, 0.94, 1.066666667e9, 11];
            ub = [0.171, 0.1845, 2.43e9, 7.5e7,        3e6, 45, 0.05, 0.15, 2.115, 2.4e9,      33];
            int_idx = 6; % n_turn
            names = {'Dp','Lgap','Kd','dP_cap','F_cap','n_turn','d_w','D_m','mu_ref','beta0','resFactor'};

        case 3
            % [b_mu, b_beta, rho_ref, alpha_rho, T_ref_C, antoine_A, antoine_B, antoine_C]
            lb = [-0.0165, -0.006, 750, 5e-4, 10, 4.5, 1200, -120];
            ub = [-0.007333333, -0.002, 950, 1.1e-3, 40, 6.5, 2200,  -20];
            int_idx = [];
            names = {'b_mu','b_beta','rho_ref','alpha_rho','T_ref_C','antoine_A','antoine_B','antoine_C'};

        case 4
            % [hA_os, hA_o_env, hA_s_env, cp_oil, cp_steel, T_env_C, T0_C, Ts0_C]
            lb = [533.3333333, 20, 20, 1500, 450, 0, 10, 10];
            ub = [1200,        200, 200, 2200, 650, 40, 40, 40];
            int_idx = [];
            names = {'hA_os','hA_o_env','hA_s_env','cp_oil','cp_steel','T_env_C','T0_C','Ts0_C'};

        case 5
            % Set 5 â€” PF + guard: [PF_tau, PF_gain, PF_t_on, p_amb, mu_min_phys, softmin_eps]
            lb = [0.30, 0.5, 0.0, 8.0e4, 0.10, 1.0e3];
            ub = [1.50, 2.5, 5.0, 1.2e5, 1.00, 1.0e5];
            int_idx = [];
            names = {'PF_tau','PF_gain','PF_t_on','p_amb','mu_min_phys','softmin_eps'};

        otherwise
            error('ga_get_bounds: Desteklenmeyen set_id=%s', num2str(set_id));
    end
end

