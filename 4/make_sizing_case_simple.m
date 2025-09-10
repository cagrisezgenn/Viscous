function [sizing, P_sized, S_worst] = make_sizing_case_simple(scaled, params, gainsPF, opts)
% A simplified, robust sizing bridge that avoids heavy reporting/printing.
% Inputs:
%   scaled, params, gainsPF (g_lo,g_mid,g_hi, PF_tau, PF_gain), opts with fields:
%     dp_allow_frac, alpha_lam, n_orf_set, dmm_step, Cd_init, mu_nom

    if nargin < 4 || isempty(opts), opts = struct(); end

    % Defaults
    dp_allow_frac = getf(opts,'dp_allow_frac',0.60);
    alpha_lam     = getf(opts,'alpha_lam',0.15);
    n_orf_set     = getf(opts,'n_orf_set',[5 6]);
    dmm_step      = getf(opts,'dmm_step',0.05);
    Cd_init       = getf(opts,'Cd_init', getf(getf(params,'orf',struct()),'CdInf',0.70));
    mu_nom        = getf(opts,'mu_nom', getf(params,'mu_ref',0.9));

    % GA lower bounds and retry controls
    d_o_lb     = 2.8e-3;                          % [m] 2.8 mm lower bound
    n_orf_lb   = max(5, min(n_orf_set));           % minimum orifice count
    Ao_lb      = n_orf_lb * (pi * d_o_lb^2 / 4);   % corresponding area
    dp_frac_min  = 0.20;                           % minimum dp_allow_frac
    dp_frac_step = 0.05;                           % step for retries

    % Lock PF and toggle gains into P
    P = params;
    nStories = size(P.M,1)-1;
    tg = ones(nStories,1) * getf(gainsPF,'g_mid',1);
    loN = min(3,nStories); if loN>0, tg(1:loN) = getf(gainsPF,'g_lo',1); end
    hiN = min(2,nStories); if hiN>0, tg(end-hiN+1:end) = getf(gainsPF,'g_hi',1); end
    P.toggle_gain = tg;
    if isfield(P,'cfg') && isfield(P.cfg,'PF')
        P.cfg.PF.tau  = getf(gainsPF,'PF_tau', getf(P.cfg.PF,'tau',1.0));
        P.cfg.PF.gain = getf(gainsPF,'PF_gain',getf(P.cfg.PF,'gain',0.85));
    end

    % Evaluate with batch windowed
    O = struct('do_export',false,'quiet',true,'thermal_reset','each','order','natural', ...
               'use_orifice',true,'use_thermal',true);
    O.mu_factors = [0.70 1.00 1.30]; O.mu_weights = [0.3 0.5 0.2];
    O.thr = Utils.default_qc_thresholds(struct());
    S_worst = run_batch_windowed(scaled, P, O);

    % Extract worst Q95 and pressure target cap
    [Q95_worst, ~, dp_cap] = find_Q95_worst_local(S_worst, O, dp_allow_frac);

    % Select orifice size and laminar length with GA bound checks
    rho = P.rho;  % density
    dp_frac = dp_allow_frac;   % working copy
    n_set   = n_orf_set;
    success = false;
    for retry = 1:10
        dp_kv_target = dp_frac * dp_cap;
        best = select_orifice_local(Q95_worst, dp_kv_target, Cd_init, n_set, dmm_step, alpha_lam, mu_nom, rho);
        best = update_Re_Cd_local(best, Q95_worst, dp_kv_target, mu_nom, P, dmm_step, rho);
        if best.d_o >= d_o_lb && best.A_o >= Ao_lb
            success = true;
            break;
        end
        if any(n_set > n_orf_lb) && best.d_o < d_o_lb
            n_set = n_orf_lb;             % try fewer orifices
            continue;
        elseif dp_frac - dp_frac_step >= dp_frac_min
            dp_frac = dp_frac - dp_frac_step;   % allow less dp, enlarge Ao
            continue;
        else
            break;
        end
    end
    if ~success
        warning('make_sizing_case_simple:bounds', ...
            'Unable to meet GA lower bounds (d_o≥%.1f mm, A_o≥%.2g m^2); enforcing minima.', ...
            d_o_lb*1e3, Ao_lb);
        best.d_o = max(best.d_o, d_o_lb);
        best.n_orf = max(best.n_orf, n_orf_lb);
        best.A_o = best.n_orf * (pi*best.d_o^2/4);
        best.Lori = alpha_lam * dp_kv_target * (pi * best.d_o^4) * best.n_orf / max(128*mu_nom*abs(Q95_worst),eps);
    end
    dp_kv_target = dp_frac * dp_cap;  % final dp target after retries

    % Build sized params
    P_sized = P;
    P_sized.orf = getf(P_sized,'orf',struct());
    P_sized.orf.d_o = best.d_o;
    P_sized.n_orf   = best.n_orf;
    P_sized.A_o     = best.A_o;
    P_sized.Lori    = best.Lori;
    P_sized.Qcap_big= max(getf(getf(P_sized,'orf',struct()),'CdInf',0.8) * P_sized.A_o, 1e-9) * sqrt(2 * 1.0e9 / P_sized.rho);
    P_sized.dp_allow_frac = dp_frac;

    % Sizing report
    sizing = struct();
    sizing.Q95_worst    = Q95_worst;
    sizing.dp_cap       = dp_cap;
    sizing.dp_kv_target = dp_kv_target;
    sizing.n_orf        = best.n_orf;
    sizing.d_o          = best.d_o;
    sizing.A_o          = best.A_o;
    sizing.Cd           = best.Cd;
    sizing.L_orif       = best.Lori;
    sizing.Qcap_big     = P_sized.Qcap_big;
    sizing.alpha_lam    = alpha_lam;
    sizing.dp_allow_frac= dp_frac;
    sizing.meets_bounds = success;
    sizing.notes        = 'Sizing (simple bridge)';
end

function [Q95_worst, dp_kv_target, dp_cap] = find_Q95_worst_local(S_worst, O, dp_allow_frac)
    Q95_worst = NaN;
    try
        vars = S_worst.table.Properties.VariableNames;
        if ismember('Q_q95_worst', vars)
            v = S_worst.table.Q_q95_worst;
            Q95_worst = max(v(:));
        elseif ismember('Q_q95_w', vars)
            v = S_worst.table.Q_q95_w;
            Q95_worst = max(v(:));
        end
    catch
    end
    dp_cap = Utils.getfield_default(O.thr,'dP95_max',1.0e9);
    dp_kv_target = dp_allow_frac * dp_cap;
end

function best = select_orifice_local(Q95_worst, dp_kv_target, Cd_init, n_orf_set, dmm_step, alpha_lam, mu_nom, rho)
    Ao_req = abs(Q95_worst) / max(Cd_init,eps) * sqrt(rho / max(2*dp_kv_target,eps));
    dgrid = @(dmm) (dmm_step * round(dmm./max(dmm_step,eps)));
    best = struct('n_orf',NaN,'d_o',NaN,'A_o',NaN,'Cd',Cd_init,'Lori',NaN);
    best_err = inf;
    for n_orf = n_orf_set(:)'
        d_o = sqrt(4*Ao_req/(pi*n_orf));
        dmm = dgrid(1000*d_o); d_o_q = dmm/1000;
        A_o = n_orf * (pi*d_o_q^2/4);
        err = abs(A_o - Ao_req)/max(Ao_req,eps);
        if err < best_err
            best_err = err; best.n_orf = n_orf; best.d_o = d_o_q; best.A_o = A_o;
        end
    end
    L_orif = alpha_lam * dp_kv_target * (pi * best.d_o^4) * best.n_orf / max(128*mu_nom*abs(Q95_worst),eps);
    best.Lori = L_orif;
end

function best = update_Re_Cd_local(best, Q95_worst, dp_kv_target, mu_nom, P, dmm_step, rho)
    dgrid = @(dmm) (dmm_step * round(dmm./max(dmm_step,eps)));
    Cd = best.Cd;
    Q_hole = Q95_worst / max(best.n_orf,1);
    Re = 4*rho*abs(Q_hole) / max(pi*mu_nom*best.d_o,eps);
    try
        Cd = P.orf.CdInf - (P.orf.CdInf - P.orf.Cd0) / (1 + (Re/max(P.orf.Rec,eps))^P.orf.p_exp);
    catch
    end
    Ao_req = abs(Q95_worst) / max(Cd,eps) * sqrt(rho / max(2*dp_kv_target,eps));
    d_o = sqrt(4*Ao_req/(pi*best.n_orf));
    dmm = dgrid(1000*d_o); best.d_o = dmm/1000; best.A_o = best.n_orf*(pi*best.d_o^2/4);
    best.Cd = Cd;
end

function v = getf(s, f, def)
    if nargin < 3, def = []; end
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = def;
    end
end

