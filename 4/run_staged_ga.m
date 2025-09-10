function out = run_staged_ga()
% RUN_STAGED_GA: 1) toggle kalibrasyonu, 2) deterministik boyutlandÄ±rma,
% 3) tohumluyup kÄ±sa GA ile rafine.

    %% Setup & data
    try, setup; catch, end
    parametreler;
    L = struct('IM_mode','band','band_fac',[0.8 1.2],'s_bounds',[0.5 2.0]);
    [~, scaled] = load_ground_motions(T1, L);

    params = struct('M',M,'C0',C0,'K',K,'k_sd',k_sd,'c_lam0',c_lam0,'Lori',Lori, ...
        'orf',orf,'rho',rho,'Ap',Ap,'A_o',A_o,'Qcap_big',Qcap_big,'mu_ref',mu_ref, ...
        'thermal',thermal,'T0_C',T0_C,'T_ref_C',T_ref_C,'b_mu',b_mu, ...
        'c_lam_min',c_lam_min,'c_lam_cap',c_lam_cap,'Lgap',Lgap, ...
        'cp_oil',cp_oil,'cp_steel',cp_steel,'steel_to_oil_mass_ratio',steel_to_oil_mass_ratio, ...
        'n_dampers_per_story',n_dampers_per_story,'toggle_gain',toggle_gain,'story_mask',story_mask, ...
        'resFactor',resFactor,'cfg',cfg,'n_orf',n_orf,'Dp',Dp,'d_w',d_w,'D_m',D_m,'n_turn',n_turn, ...
        'Kd',Kd,'Ebody',Ebody,'Gsh',Gsh,'story_height',story_height);

    %% 1) Toggle kalibrasyonu
    T = run_toggle_calibration(scaled, params);
    gstar = T.g_star;
    fprintf('toggle_gain* = %.3f\n', gstar);

    %% 2) Deterministik boyutlandÄ±rma (kÃ¶prÃ¼)
    gainsPF = struct('g_lo',gstar,'g_mid',gstar,'g_hi',gstar, ...
                     'PF_tau', Utils.getfield_default(params.cfg.PF,'tau',1.0), ...
                     'PF_gain',Utils.getfield_default(params.cfg.PF,'gain',0.85));
    Sopt = struct('dp_allow_frac',0.45,'alpha_lam',0.15,'n_orf_set',[3 4 5], 'dmm_step',0.1,'Cd_init',0.75,'report_toggle',false);
    [sizing, P_sized, S_worst] = make_sizing_case_simple(scaled, params, gainsPF, Sopt); %#ok<ASGLU>

    %% 3) GA iÃ§in tohum ve kÄ±sa rafine
    step_vec = [0.1 NaN 0.01 0.02 0.01 0.01 0.05 1 25 1 0.5 5 NaN 0.05];
    lb = [2.80, 5, 0.95, 0.78, 0.50, 0.75, 0.80, 60, 200, 100, 8, 60, 6, 0.60];
    ub = [3.60, 8, 1.10, 0.90, 0.70, 0.95, 1.40, 140, 800, 160, 16, 100, 12, 1.20];

    % Tohum parametreleri GA sÄ±nÄ±rlarÄ± iÃ§inde mi kontrol et
    x0 = encode_x_from_params(P_sized);
    out_of_range = (x0 < lb) | (x0 > ub);
    if any(out_of_range)
        lb(out_of_range) = min(lb(out_of_range), x0(out_of_range));
        ub(out_of_range) = max(ub(out_of_range), x0(out_of_range));
        fprintf('Deterministic design outside GA bounds; expanding bounds.\n');
    end

    % KomÅŸuluk Ã¶rnekleri (nicemli, sÄ±nÄ±rda kÄ±rpÄ±lmÄ±ÅŸ)
    % Komşuluk örnekleri (nicemli, sınırda kırpılmış) — yerel fonksiyon qclamp_local kullanılır
    Nseed = 40; P0 = zeros(Nseed, numel(x0));
    P0(1,:) = qclamp_local(x0, lb, ub);
    rng(42);
    for i=2:Nseed
        xr = x0 + (randn(size(x0)).* [0.1 0 0.01 0.02 0.01 0.01 0.05 2 50 2 0.5 5 0 0.05]);
        xr(2) = round(x0(2) + randi([-1 1]));
        xr(13)= round(x0(13) + randi([-1 1]));
        P0(i,:) = qclamp_local(xr, lb, ub);
    end

    optsEval = struct('mu_factors',[0.70 1.00 1.30],'mu_weights',[0.3 0.5 0.2], ...
                      'penalty_scale',5,'penalty_power',1.0, ...
                      'penalty_weights',struct('dP',1,'Qcap',1,'cav',3,'T',1,'mu',0.5));
    optsGA = struct('PopulationSize', Nseed, 'MaxGenerations', 40, 'CrossoverFraction',0.8, ...
                    'ParetoFraction',0.5, 'UseParallel', true, 'InitialPopulationMatrix', P0);

    [X,F,gaout] = run_ga_driver(scaled, params, optsEval, optsGA); %#ok<ASGLU>

    out = struct('g_star',gstar,'x0',x0);
end

function x = encode_x_from_params(P)
    x = nan(1,14);
    x(1)  = 1e3*P.orf.d_o;
    x(2)  = P.n_orf;
    x(3)  = Utils.getfield_default(P.cfg.PF,'tau',1.0);
    x(4)  = Utils.getfield_default(P.cfg.PF,'gain',0.85);
    x(5)  = P.orf.Cd0;
    x(6)  = P.orf.CdInf;
    x(7)  = P.orf.p_exp;
    x(8)  = 1e3*P.Lori;
    x(9)  = Utils.getfield_default(P.thermal,'hA_W_perK',450);
    x(10) = 1e3*P.Dp;
    x(11) = 1e3*P.d_w;
    x(12) = 1e3*P.D_m;
    x(13) = P.n_turn;
    x(14) = P.mu_ref;
end



function xq = qclamp_local(x, lb, ub)
    xq = x;
    xq(1) = Utils.quantize_step(xq(1),0.05);
    xq(3) = Utils.quantize_step(xq(3),0.01);
    xq(4) = Utils.quantize_step(xq(4),0.02);
    xq(5) = Utils.quantize_step(xq(5),0.01);
    xq(6) = Utils.quantize_step(xq(6),0.01);
    xq(7) = Utils.quantize_step(xq(7),0.05);
    xq(8) = Utils.quantize_step(xq(8),1);
    xq(9) = Utils.quantize_step(xq(9),25);
    xq(10)= Utils.quantize_step(xq(10),1);
    xq(11)= Utils.quantize_step(xq(11),0.5);
    xq(12)= Utils.quantize_step(xq(12),5);
    xq(13)= round(xq(13));
    xq(14)= Utils.quantize_step(xq(14),0.05);
    xq = max(lb, min(ub, xq));
    xq(13)=round(xq(13));
end
