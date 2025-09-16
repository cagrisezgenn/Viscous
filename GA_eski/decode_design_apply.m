function [geom, sh, orf, hyd, therm, num, ga, env] = decode_design_apply(ga, geom, sh, orf, hyd, therm, num, env)
%DECODE_DESIGN_APPLY  Decode GA design vector and apply to model structures.
%   The design variables are described by a table returned from
%   GET_DESIGN_TABLE.  Each entry provides the variable name, bounds,
%   integer flag, target structure/field and an optional post-processing
%   rule.  This function clamps the GA vector to bounds, applies it to the
%   target structures and finally executes any post rules (such as
%   cross-variable checks).

    if nargin < 8
        env = struct();
    end
    if ~isfield(env, 'LOG') || isempty(env.LOG)
        env.LOG = struct();
    end

    % set_id güvenliği
    set_id = ga.design_set;
    if isstring(set_id) || ischar(set_id), set_id = str2double(set_id); end
    if isempty(set_id) || isnan(set_id), set_id = 1; end
    ga.design_set = double(set_id);

    % Tasarım değişken tablolarını al
    vars   = get_design_table(ga.design_set);
    lb     = [vars.lb];
    ub     = [vars.ub];
    names  = {vars.name};
    int_idx = find([vars.int]);

    % GA kapalıysa: sadece bounds/meta döndür
    if ~isfield(ga,'enable') || ~ga.enable || isempty(ga.x)
        ga.lb=lb; ga.ub=ub; ga.int_idx=int_idx; ga.names=names; ga.x_use=[];
        return;
    end

    x = ga.x(:);
    if numel(x)~=numel(lb)
        error('ga.x uzunluğu set-%d için %d olmalı.', ga.design_set, numel(lb));
    end

    % Kenet + tamsayı
    x = max(lb(:), min(ub(:), x(:)));
    if ~isempty(int_idx)
        x(int_idx) = round(x(int_idx));
        x = max(lb(:), min(ub(:), x(:)));
    end

    % Değerleri ilgili strukturlara uygula
    for ii = 1:numel(vars)
        val = x(ii);
        switch vars(ii).target
            case 'geom'
                geom.(vars(ii).field) = val;
            case 'sh'
                sh.(vars(ii).field) = val;
            case 'orf'
                orf.(vars(ii).field) = val;
            case 'hyd'
                hyd.(vars(ii).field) = val;
            case 'therm'
                therm.(vars(ii).field) = val;
            case 'num'
                num.(vars(ii).field) = val;
            case 'ga'
                ga.(vars(ii).field) = val;
            case 'env'
                env.(vars(ii).field) = val;
            otherwise
                error('decode_design_apply: bilinmeyen hedef %s', vars(ii).target);
        end
    end

    % Son işlem kurallarını uygula
    for ii = 1:numel(vars)
        if ~isempty(vars(ii).post)
            vars(ii).post();
        end
    end

    % Hat ataletini geometri ve yoğunluktan türeterek hesapla
    Ao = pi*geom.d_o^2/4;
    hyd.Lh = therm.rho_ref * geom.Lori / ( (orf.n_orf * Ao)^2 );

    ga.lb=lb; ga.ub=ub; ga.int_idx=int_idx; ga.names=names; ga.x_use=x;
    try
        if isfield(env.LOG,'verbose_decode') && env.LOG.verbose_decode
            fprintf('GA decode: set=%s | x_use = [%s]\n', num2str(ga.design_set), strjoin(compose('%.3g',x.'),', '));
        end
    catch
    end

    % ==== Nested functions ==================================================

    function vars = get_design_table(set_id)
        switch set_id
            case 1
                % Set 1: [Cd0, CdInf, Rec, p_exp, cav_sf, K_leak, Vmin_fac, d_o, Lori]
                vars(1) = make_var('Cd0',      0.30, 0.90, 0, 'orf',  'Cd0',      []);
                vars(2) = make_var('CdInf',    0.60, 1.35, 0, 'orf',  'CdInf',    @post_CdInf);
                vars(3) = make_var('Rec',   2533.333333, 5700.000000, 0, 'orf',  'Rec',      []);
                vars(4) = make_var('p_exp',    0.55, 1.65, 0, 'orf',  'p_exp',    []);
                vars(5) = make_var('cav_sf',   0.525,1.575,0, 'orf',  'cav_sf',   []);
                vars(6) = make_var('K_leak',   1.0e-9,6.0e-9,0,'hyd','K_leak',   []);
                vars(7) = make_var('Vmin_fac', 0.93, 0.98, 0,'hyd','Vmin_fac', []);
                vars(8) = make_var('d_o',      1.0e-4,5.0e-3,0,'geom','d_o', []);
                vars(9) = make_var('Lori',     1.0e-4,1.0e-2,0,'geom','Lori', []);
            case 2
                % Set 2: [Dp, Lgap, Kd, dP_cap, F_cap, n_turn, d_w, D_m, mu_ref, beta0, resFactor]
                vars(1) = make_var('Dp',       0.076, 0.171, 0,'geom','Dp',       []);
                vars(2) = make_var('Lgap',     0.0615,0.1845,0,'geom','Lgap',     []);
                vars(3) = make_var('Kd',       1.08e9,2.43e9,0,'geom','Kd',       []);
                vars(4) = make_var('dP_cap',   3.333333333e7,7.5e7,0,'num','dP_cap',[]);
                vars(5) = make_var('F_cap',    1e6,   3e6,   0,'num','F_cap',    []);
                vars(6) = make_var('n_turn',   15,    45,    1,'sh', 'n_turn',    []);
                vars(7) = make_var('d_w',      0.01, 0.05, 0,'sh','d_w', []);
                vars(8) = make_var('D_m',      0.05, 0.15, 0,'sh','D_m', []);
                vars(9) = make_var('mu_ref',   0.94,  2.115, 0,'therm','mu_ref', []);
                vars(10) = make_var('beta0',   1.066666667e9,2.4e9,0,'therm','beta0', []);
                vars(11) = make_var('resFactor',11,    33,    0,'therm','resFactor', []);
            case 3
                % Set 3: [b_mu, b_beta, rho_ref, alpha_rho, T_ref_C, antoine_A, antoine_B, antoine_C]
                vars(1) = make_var('b_mu',     -0.0165, -0.007333333,0,'therm','b_mu', []);
                vars(2) = make_var('b_beta',   -0.006,  -0.002,      0,'therm','b_beta', []);
                vars(3) = make_var('rho_ref',  750,    950,          0,'therm','rho_ref', []);
                vars(4) = make_var('alpha_rho',5e-4,  1.1e-3,        0,'therm','alpha_rho', []);
                vars(5) = make_var('T_ref_C',  10,     40,           0,'therm','T_ref_C', []);
                vars(6) = make_var('antoine_A',4.5,    6.5,          0,'therm','antoine_A', []);
                vars(7) = make_var('antoine_B',1200,   2200,         0,'therm','antoine_B', []);
                vars(8) = make_var('antoine_C',-120,   -20,          0,'therm','antoine_C', []);
            case 4
                % Set 4: [hA_os, hA_o_env, hA_s_env, cp_oil, cp_steel, T_env_C, T0_C, Ts0_C]
                vars(1) = make_var('hA_os',    533.3333333, 1200, 0,'therm','hA_os', []);
                vars(2) = make_var('hA_o_env', 20, 200, 0,'therm','hA_o_env', []);
                vars(3) = make_var('hA_s_env', 20, 200, 0,'therm','hA_s_env', []);
                vars(4) = make_var('cp_oil',   1500, 2200, 0,'therm','cp_oil', []);
                vars(5) = make_var('cp_steel', 450, 650, 0,'therm','cp_steel', []);
                vars(6) = make_var('T_env_C',  0,   40,  0,'therm','T_env_C', []);
                vars(7) = make_var('T0_C',     10,  40,  0,'therm','T0_C', []);
                vars(8) = make_var('Ts0_C',    10,  40,  0,'therm','Ts0_C', []);
            case 5
                % Set 5: [PF_tau, PF_gain, PF_t_on, p_amb, mu_min_phys, softmin_eps]
                vars(1) = make_var('PF_tau',    0.30, 1.50, 0,'ga',  'PF_tau',    []);
                vars(2) = make_var('PF_gain',   0.5,  2.5,  0,'ga',  'PF_gain',   []);
                vars(3) = make_var('PF_t_on',   0.0,  5.0,  0,'ga',  'PF_t_on',   []);
                vars(4) = make_var('p_amb',     8.0e4,1.2e5,0,'orf', 'p_amb',     []);
                vars(5) = make_var('mu_min_phys',0.10,1.00,0,'num', 'mu_min_phys', []);
                vars(6) = make_var('softmin_eps',1.0e3,1.0e5,0,'num','softmin_eps', []);
            otherwise
                error('decode_design_apply: Desteklenmeyen design_set=%s', num2str(set_id));
        end
    end

    function v = make_var(name, lb, ub, isint, target, field, post)
        v = struct('name',name,'lb',lb,'ub',ub,'int',logical(isint), ...
                   'target',target,'field',field,'post',post);
    end

    function post_CdInf()
        if orf.CdInf < orf.Cd0
            orig = orf.CdInf;
            orf.CdInf = max(orf.Cd0, orf.CdInf);
            error('decode_design_apply:CdInfBound', ...
                  'CdInf (%.3g) < Cd0 (%.3g). Ayarlandı: %.3g.', ...
                  orig, orf.Cd0, orf.CdInf);
        end
    end

end

