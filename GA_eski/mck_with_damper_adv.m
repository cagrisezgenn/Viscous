function [x,a,diag,varargout] = mck_with_damper_adv(t,ag,M,C,K,k_sd,geom,orf,hyd,therm,num,cfg)
    nd = 1; if isfield(hyd,'n_parallel'), nd = max(1, round(hyd.n_parallel)); end
    Ao1 = max(orf.n_orf * (pi*geom.d_o^2/4), 1e-12);
    Ao  = nd * Ao1;              % toplam orifis alanÄ±
    Ap1 = geom.Ap;
    Ap  = nd * Ap1;              % toplam piston alanÄ±

    n = size(M,1); r = ones(n,1); Ns = n-1;
    z0 = zeros(2*n + 2*Ns + Ns + 2, 1);
    z0(2*n + (1:Ns)) = orf.p_amb;  z0(2*n + Ns + (1:Ns)) = orf.p_amb;
    z0(end-1)=therm.T0_C; z0(end)=therm.Ts0_C;

    % --- Ã–lÃ§ekli toleranslar (basÄ±nca ve akÄ±ÅŸa gÃ¶re) ---
    rho0  = therm.rho_ref;
    dPcap = max(num.dP_cap, 1e5);
    Qcap_est = Ao * sqrt( 2*dPcap / max(rho0,100) );   % Ao = nd*Ao1

    AbsX = 1e-4;  AbsV = 1e-3;  AbsP = max(5e3, 1e-5 * dPcap);
    AbsQ = max(1e-5, 1e-3 * Qcap_est);  AbsT = 5e-3;
    AbsTol_default = [AbsX*ones(n,1); AbsV*ones(n,1); AbsP*ones(2*Ns,1); AbsQ*ones(Ns,1); AbsT*ones(2,1)];

    % ----- Solver tolerans/step seÃ§enekleri (cfg.solver Ã¼zerinden) -----
    dt  = median(diff(t));
    % Preset â†’ taban deÄŸerler
    preset = 'balanced';
    if isfield(cfg,'solver') && isfield(cfg.solver,'preset') && ~isempty(cfg.solver.preset)
        preset = lower(string(cfg.solver.preset));
    end
    switch char(preset)
        case 'explore'
            % Daha gevÅŸek tolerans + daha temkinli adÄ±mlar
            RelTol_p = 1.0e-1; MaxStep_p = max(12*dt, 3e-3);  InitialStep_p = max(0.25*dt, 1.0e-3);
        case 'tight'
            % SÄ±kÄ± tolerans + kÃ¼Ã§Ã¼k adÄ±mlar
            RelTol_p = 3.0e-2; MaxStep_p = max(8*dt,  2e-3);  InitialStep_p = max(0.20*dt, 1.0e-3);
        otherwise % 'balanced'
            % Dengeli tolerans + gÃ¼venli adÄ±mlar
            RelTol_p = 7.0e-2; MaxStep_p = max(12*dt, 3e-3);  InitialStep_p = max(0.25*dt, 1.5e-3);
    end
    % Override'lar
    if isfield(cfg,'solver') && isstruct(cfg.solver)
        if isfield(cfg.solver,'RelTol') && isfinite(cfg.solver.RelTol)
            RelTol_p = cfg.solver.RelTol;
        end
        if isfield(cfg.solver,'MaxStep') && isfinite(cfg.solver.MaxStep) && cfg.solver.MaxStep>0
            MaxStep_p = cfg.solver.MaxStep;
        end
        if isfield(cfg.solver,'InitialStep') && isfinite(cfg.solver.InitialStep) && cfg.solver.InitialStep>0
            InitialStep_p = cfg.solver.InitialStep;
        end
    end
    % AbsTol: varsayÄ±lan vektÃ¶r + Ã¶lÃ§ek/override
    AbsTol = AbsTol_default;
    if isfield(cfg,'solver') && isstruct(cfg.solver)
        if isfield(cfg.solver,'AbsTol') && ~isempty(cfg.solver.AbsTol)
            at = cfg.solver.AbsTol;
            if isnumeric(at)
                if isscalar(at) && isfinite(at) && at>0
                    AbsTol = AbsTol * at;  % skaler: Ã¶lÃ§ek faktÃ¶rÃ¼
                elseif numel(at)==numel(AbsTol_default) && all(isfinite(at(:)))
                    AbsTol = at(:);
                end
            end
        end
        if isfield(cfg.solver,'AbsTolScale') && isfinite(cfg.solver.AbsTolScale) && cfg.solver.AbsTolScale>0
            AbsTol = AbsTol * cfg.solver.AbsTolScale;
        end
    end

    opts = odeset('RelTol',RelTol_p,'AbsTol',AbsTol,'JPattern',local_JPattern(n), ...
                  'MaxStep',MaxStep_p,'InitialStep',InitialStep_p);

    agf = griddedInterpolant(t, ag, 'pchip', 'nearest');

    odef=@(tt,z) rhs(tt,z);
    sol=ode23tb(odef,[t(1) t(end)],z0,opts);

    t_end_sol=sol.x(end); t_use=t(t<=t_end_sol+1e-12); Z_use=deval(sol,t_use).';
    if t_end_sol < t(end)-1e-9
        % fallback: biraz daha gevÅŸek RelTol kullan (preset tabanÄ±nÄ±n bir kademe Ã¼stÃ¼)
        opts2 = odeset(opts, 'RelTol', max(RelTol_p*1.2, RelTol_p), 'MaxOrder', 2);
        try
            sol2=ode15s(odef,[t_end_sol t(end)],sol.y(:,end),opts2);
            t_use2=t(t>t_end_sol); Z_use2=deval(sol2,t_use2).'; Z=[Z_use;Z_use2];
            warning('ode23tb early stop at t=%.3f â†’ continued to t=%.3f.',t_end_sol,t(end));
        catch
            Z=[Z_use; nan(numel(t)-numel(t_use), size(Z_use,2))];
            warning('Both solvers stopped early at t=%.3f s.',t_end_sol);
        end
    else
        Z=Z_use;
    end

    x = Z(:,1:n); v = Z(:,n+1:2*n);
    Fdev = dev_force_from_story(t, x, Z, k_sd, geom, hyd, therm, num, orf, cfg);
    a = ( -(M\(C*v.' + K*x.' + Fdev)).' - ag.*r.' );

    [drift, F_story, dP_orf_t, T_oil, T_steel, mu_t, E_cum, ...
     cav_frac_t, dP_q50, dP_q95, Q_abs_med, Q_abs_p95, ...
     cav_margin_t, cav_margin_min] = ...
        diagnostics_time_series(t, Z, n, Ns, k_sd, geom, orf, hyd, therm, num, cfg);

    diag = struct('drift',drift,'F_story',F_story,'dP_orf_time',dP_orf_t, ...
                  'dP_orf_time_max',max(dP_orf_t,[],2),'T_oil',T_oil,'T_steel',T_steel, ...
                  'mu_t',mu_t,'E_cum',E_cum,'cav_frac_t',cav_frac_t,'dP_q50',dP_q50,'dP_q95',dP_q95, ...
                  'Q_abs_med',Q_abs_med,'Q_abs_p95',Q_abs_p95, ...
                  'cav_margin_t',cav_margin_t,'cav_margin_min',cav_margin_min);

    if nargout>=4, varargout{1} = v; end

    % -------------------- iÃ§ fonksiyonlar --------------------
    function dz = rhs(tt,z)
        x  = z(1:n); v = z(n+1:2*n);
        p1 = z(2*n + (1:Ns)); p2 = z(2*n + Ns + (1:Ns));
        Q  = z(2*n + 2*Ns + (1:Ns)); T_o = z(end-1); T_s = z(end);

        mu_raw = therm.mu_ref * exp(therm.b_mu*(T_o - therm.T_ref_C));
        if cfg.use_thermal && cfg.on.mu_floor, mu=max(num.mu_min_phys,mu_raw); else, mu=cfg.use_thermal*mu_raw + (~cfg.use_thermal)*therm.mu_ref; end
        if cfg.use_thermal
            rho=max(100,therm.rho_ref/(1+therm.alpha_rho*(T_o-therm.T_ref_C)));
            beta=max(1e8,therm.beta0*exp(therm.b_beta*(T_o-therm.T_ref_C)));
            p_vap=p_vap_Antoine(T_o,therm,orf);
        else
            rho=therm.rho_ref; beta=therm.beta0; p_vap=p_vap_Antoine(therm.T_ref_C,therm,orf);
        end

        drift = x(2:end) - x(1:end-1);
        dvel  = v(2:end) - v(1:end-1);

        V1 = nd*hyd.V0 + Ap*drift;
        V2 = nd*hyd.V0 - Ap*drift;
        Vmin = nd * hyd.Vmin_fac * hyd.V0; V1=max(V1,Vmin); V2=max(V2,Vmin);

        % --- Orifis hidrolik kayÄ±plarÄ± ---
        if cfg.use_orifice
            if cfg.on.CdRe
                Re = rho .* abs(Q) .* geom.d_o ./ max(Ao * mu, 1e-9);
                Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
                Cd = max(min(Cd, 1.2), 0.2);
            else
                Cd = orf.CdInf;
            end
            if cfg.on.Rkv, RQ = rho ./ max(2 * (Cd .* Ao).^2, 1e-12); else, RQ = 0 * Q; end

            Qcap = getfield_default(num,'Qcap_big', ...
                0.4 * ( max(max(orf.CdInf,orf.Cd0)*Ao, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
            if cfg.on.Qsat, Q_h = Qcap * tanh(Q ./ max(Qcap, 1e-9)); else, Q_h = Q; end

            dP_kv = RQ .* Q_h .* abs(Q_h);
            if cfg.on.Rlam
                R_lam  = (128 * mu * geom.Lori / (pi * geom.d_o^4)) / max(1, nd * orf.n_orf);
                dP_lam = R_lam .* Q;
            else
                dP_lam = 0 * Q;
            end
            dP_h = dP_lam + dP_kv;
        else
            dP_h = 0 * Q; Q_h = 0 * Q;
        end

        % --- Kavitasyon, dP_cap, atalet, leak ---
        if cfg.on.cavitation, p2_eff = max(p2, orf.cav_sf * p_vap); else, p2_eff = p2; end
        dP_raw = p1 - p2_eff;
        if cfg.on.dP_cap && isfinite(num.dP_cap)
            dP_eff = num.dP_cap * tanh(dP_raw ./ max(num.dP_cap, 1));
        else
            dP_eff = dP_raw;
        end
        if cfg.use_orifice && cfg.on.hyd_inertia, dQ = (dP_eff - dP_h) ./ max(hyd.Lh, 1e-12); else, dQ = 0 * Q; end
        if cfg.on.leak, Q_leak = hyd.K_leak * (p1 - p2); else, Q_leak = 0 * Q; end

        % --- BasÄ±nÃ§ ODE'leri ---
        if cfg.on.pressure_ode
            dp1 = (beta ./ V1) .* ( -Q - Q_leak - Ap * dvel );
            dp2 = (beta ./ V2) .* ( +Q + Q_leak + Ap * dvel );
        else
            dp1 = 0 * p1; dp2 = 0 * p2;
        end
        % Cavitation clamp
        if cfg.on.cavitation
            m1 = (p1 <= p_vap) & ((-Q - Q_leak - Ap*dvel) < 0);
            m2 = (p2 <= p_vap) & ((+Q + Q_leak + Ap*dvel) < 0);
            dp1(m1) = 0; dp2(m2) = 0;
        end

        % --- Damper kuvveti (yay + PF) â€” RESISTIVE-ONLY KLAMP BURADA ---
        w_pf = cfg.on.pressure_force * (1 - exp(-max(tt - cfg.PF.t_on, 0) / max(cfg.PF.tau, 1e-6)));
        F_story = k_sd * drift;

        if cfg.on.pressure_force
            dp_pf = (p1 - p2_eff);                % Ns x 1
            if (isfield(cfg,'PF') && isfield(cfg.PF,'resistive_only') && cfg.PF.resistive_only) || (isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only)
s = tanh(20*dvel);   % 20 ~ 1/(0.05 m/s) benzeri yumuÅŸatma
                % Ä°stersen daha yumuÅŸak iÃ§in: s = tanh(20*dvel);
                dp_pf = s .* max(0, s .* dp_pf);  % yalnÄ±z direnÃ§li bileÅŸen
            end
            F_story = F_story + Ap * (cfg.PF.gain * w_pf .* dp_pf);  % w_pf skaler
        end

        % --- Kuvvet daÄŸÄ±tÄ±mÄ±, yapÄ± ODE ---
        F        = zeros(n,1);
        F(1)     = -F_story(1);
        if n > 2, F(2:n-1) = F_story(1:end-1) - F_story(2:end); end
        F(n)     =  F_story(end);

        dv = M \ ( -C*v - K*x - F - M*r*agf(tt) );

        % --- IsÄ± ODE'leri (toplam kayÄ±p gÃ¼cÃ¼) ---
        if cfg.use_thermal
            if cfg.use_orifice && cfg.on.Rlam, P_lam = dP_lam .* Q; else, P_lam = 0; end
            if cfg.use_orifice && cfg.on.Rkv,  P_kv  = dP_kv  .* Q; else, P_kv  = 0; end
            P_sum = sum(P_lam + P_kv); P_sum = max(P_sum, 0);

            dT_o = ( P_sum ...
                     - therm.hA_os   * (T_o - T_s) ...
                     - therm.hA_o_env* (T_o - therm.T_env_C) ) / max(therm.C_oil,   eps);
            dT_s = ( + therm.hA_os   * (T_o - T_s) ...
                     - therm.hA_s_env* (T_s - therm.T_env_C) ) / max(therm.C_steel, eps);
        else
            dT_o = 0; dT_s = 0;
        end

        dz = [ v; dv; dp1; dp2; dQ; dT_o; dT_s ];
    end

  function Jp = local_JPattern(nn)
    % Bant (seyrek) Jacobian ÅŸablonu: sadece gerÃ§ek etkileÅŸimler.
    % Bloklar: [x(1..n), v(1..n), p1(1..Ns), p2(1..Ns), Q(1..Ns), T_o, T_s]
    % SÄ±ralama rhs(): [v; dv; dp1; dp2; dQ; dT_o; dT_s]

    % Boyutlar
    Ns   = nn - 1;
    Ntot = 2*nn + 2*Ns + Ns + 2;  % 2n (x,v) + 2Ns (p1,p2) + Ns (Q) + 2 (T_o,T_s)

    % Blok indeksleri
    ix  = 1:nn;                     % x
    iv  = nn + (1:nn);              % v
    ip1 = 2*nn + (1:Ns);            % p1
    ip2 = 2*nn + Ns + (1:Ns);       % p2
    iQ  = 2*nn + 2*Ns + (1:Ns);     % Q
    iTo = 2*nn + 3*Ns + 1;          % T_o
    iTs = 2*nn + 3*Ns + 2;          % T_s

    I = []; J = [];  % (row, col) Ã§iftleri

    % 1) x' = v  -> her i iÃ§in (row=ix(i), col=iv(i))
    I = [I, ix];
    J = [J, iv];

    % 2) v'  â€” yapÄ±sal sistem + damper basÄ±ncÄ±
    % - x komÅŸuluk (i-1,i,i+1) Ã¼zerinden
    % - v komÅŸuluk (i-1,i,i+1) Ã¼zerinden (C ve dvel)
    % - p1/p2 hikaye(i-1) ve hikaye(i) Ã¼zerinden
    % - T_o (p_vap ile PF baÄŸlamÄ±) Ã¼zerinden
    for i = 1:nn
      row = iv(i);

      % x komÅŸuluklarÄ±
      xn = i + [-1 0 1]; xn = xn(xn>=1 & xn<=nn);
      I = [I, row*ones(1, numel(xn))];
      J = [J, ix(xn)];

      % v komÅŸuluklarÄ±
      vn = i + [-1 0 1]; vn = vn(vn>=1 & vn<=nn);
      I = [I, row*ones(1, numel(vn))];
      J = [J, iv(vn)];

      % p1/p2 hikaye komÅŸuluklarÄ±
      if i==1
        ks = 1;
      elseif i==nn
        ks = Ns;
      else
        ks = [i-1, i];
      end
      I = [I, row*ones(1, numel(ks)), row*ones(1, numel(ks))];
      J = [J, ip1(ks),                     ip2(ks)                    ];

      % T_o etkisi (tek sÃ¼tun)
      I = [I, row];
      J = [J, iTo];
    end

    % 3) p1'  â€” her hikaye k, yerel x(k:k+1), v(k:k+1), p1(k), p2(k), Q(k), T_o
    for k = 1:Ns
      row = ip1(k);
      % x(k), x(k+1)
      I = [I, row, row];
      J = [J, ix(k), ix(k+1)];
      % v(k), v(k+1)
      I = [I, row, row];
      J = [J, iv(k), iv(k+1)];
      % p1(k), p2(k)
      I = [I, row, row];
      J = [J, ip1(k), ip2(k)];
      % Q(k)
      I = [I, row];
      J = [J, iQ(k)];
      % T_o
      I = [I, row];
      J = [J, iTo];
    end

    % 4) p2'  â€” her hikaye k, yerel x(k:k+1), v(k:k+1), p1(k), p2(k), Q(k), T_o
    for k = 1:Ns
      row = ip2(k);
      % x(k), x(k+1)
      I = [I, row, row];
      J = [J, ix(k), ix(k+1)];
      % v(k), v(k+1)
      I = [I, row, row];
      J = [J, iv(k), iv(k+1)];
      % p1(k), p2(k)
      I = [I, row, row];
      J = [J, ip1(k), ip2(k)];
      % Q(k)
      I = [I, row];
      J = [J, iQ(k)];
      % T_o
      I = [I, row];
      J = [J, iTo];
    end

    % 5) Q'   â€” her hikaye k, p1(k), p2(k), Q(k), T_o
    for k = 1:Ns
      row = iQ(k);
      I = [I, row, row, row, row];
      J = [J, ip1(k), ip2(k), iQ(k), iTo];
    end

    % 6) T_o' â€” tÃ¼m Q'lar, T_o, T_s
    I = [I, iTo*ones(1,Ns), iTo, iTo];
    J = [J, iQ,                iTo, iTs];

    % 7) T_s' â€” T_o, T_s
    I = [I, iTs, iTs];
    J = [J, iTo, iTs];

    % (row,col) -> seyrek 1'ler
    Jp = sparse(I, J, 1, Ntot, Ntot);
  end

 end

function [drift, F_story, dP_orf_t, T_oil, T_steel, mu_t, E_cum, ...
          cav_frac_t, dP_q50, dP_q95, Q_abs_med, Q_abs_p95, ...
          cav_margin_t, cav_margin_min] = ...
          diagnostics_time_series(t, Z, n, Ns, k_sd, geom, orf, hyd, therm, num, cfg)

    nd  = 1; if isfield(hyd,'n_parallel'), nd = hyd.n_parallel; end
    Ao1 = orf.n_orf * (pi*geom.d_o^2/4);
    Ao  = nd * Ao1;
    Ap  = nd * geom.Ap;

    X  = Z(:,1:n);
    V  = Z(:,n+1:2*n);
    p1 = Z(:,2*n + (1:Ns));
    p2 = Z(:,2*n + Ns + (1:Ns));
    Q  = Z(:,2*n + 2*Ns + (1:Ns));
    T_oil   = Z(:,end-1);
    T_steel = Z(:,end);

    drift = X(:,2:end) - X(:,1:end-1);
    dvel  = V(:,2:end) - V(:,1:end-1);      % Nt x Ns

    mu_raw  = therm.mu_ref * exp(therm.b_mu*(T_oil - therm.T_ref_C));
    if cfg.use_thermal && cfg.on.mu_floor, mu_t = max(num.mu_min_phys, mu_raw);
    else,                                   mu_t = cfg.use_thermal.*mu_raw + (~cfg.use_thermal).*therm.mu_ref; end

    rho_t   = max(100, therm.rho_ref ./ (1 + therm.alpha_rho*(T_oil - therm.T_ref_C)));
    p_vap_t = p_vap_Antoine(T_oil, therm, orf);

    cav_margin_t  = min(p2 - p_vap_t, [], 2);
    cav_margin_min = min(cav_margin_t, [], 'omitnan');

    if cfg.on.cavitation, p2_eff = max(p2, orf.cav_sf * p_vap_t); else, p2_eff = p2; end

    % ---- Hidrolik kayÄ±plar (dP_h) + dP_orf_t ----
    if cfg.use_orifice
        if cfg.on.CdRe
            Re = rho_t .* abs(Q) .* geom.d_o ./ max(Ao .* mu_t, 1e-9);
            Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
            Cd = max(min(Cd, 1.2), 0.2);
        else
            Cd = orf.CdInf;
        end
        if cfg.on.Rkv, RQ = rho_t ./ max(2 * (Cd .* Ao).^2, 1e-12); else, RQ = 0 * Q; end

        Qcap = getfield_default(num,'Qcap_big', ...
            0.4 * ( max(max(orf.CdInf,orf.Cd0)*Ao, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
        if cfg.on.Qsat, Q_h = Qcap * tanh(Q ./ max(Qcap, 1e-9)); else, Q_h = Q; end

        dP_kv = RQ .* Q_h .* abs(Q_h);
        if cfg.on.Rlam
            R_lam_t = (128 * mu_t .* geom.Lori ./ (pi * geom.d_o^4)) / max(1, nd * orf.n_orf);
            dP_lam  = R_lam_t .* Q;
        else
            dP_lam  = 0 * Q;
        end
        dP_h = dP_lam + dP_kv;
    else
        Q_h = 0 * Q; dP_h = 0 * Q;
    end

    dP_raw = p1 - p2_eff;
    if cfg.on.dP_cap && isfinite(num.dP_cap)
        dP_eff = num.dP_cap * tanh( dP_raw ./ max(num.dP_cap, 1) );
    else
        dP_eff = dP_raw;
    end
    epsm     = max(1e3, double(num.softmin_eps));
    dP_orf_t = 0.5 * ( dP_eff + dP_h - sqrt( (dP_eff - dP_h).^2 + epsm^2 ) );

    % ---- PF kuvveti: RESISTIVE-ONLY burada uygulanÄ±yor ----
    w_pf_vec = pf_weight(t, cfg) * cfg.PF.gain;
    if cfg.on.pressure_force
        dP_pf = dP_orf_t;                               % Nt x Ns
        if (isfield(cfg,'PF') && isfield(cfg.PF,'resistive_only') && cfg.PF.resistive_only) || (isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only)
            s = tanh(20*dvel);   % 20 ~ 1/(0.05 m/s) benzeri yumuÅŸatma
            dP_pf = s .* max(0, s .* dP_pf);
        end
        F_story = k_sd*drift + (w_pf_vec .* dP_pf) * Ap;  % Nt x Ns â†’ (.*) yayÄ±mlÄ±
    else
        F_story = k_sd*drift;
    end

    cav_mask   = p2 < p_vap_t;
    cav_frac_t = mean(cav_mask, 2);

    if cfg.use_orifice
        P_lam = dP_lam .* Q; P_kv  = dP_kv  .* Q;
        P_sum = sum(P_lam + P_kv, 2); P_sum = max(P_sum, 0);
    else
        P_sum = zeros(size(t));
    end
    E_cum = cumtrapz(t, P_sum);

    dP_q50    = prctile(dP_orf_t, 50, 2);
    dP_q95    = prctile(dP_orf_t, 95, 2);
    qQ        = prctile(abs(Q(:)), [50 95]);
    Q_abs_med = qQ(1);
    Q_abs_p95 = qQ(2);
end

function F = dev_force_from_story(t, X, Z, k_sd, geom, hyd, therm, num, orf, cfg)
    Nt=size(X,1); n=size(X,2); Ns=n-1;
    drift = X(:,2:end) - X(:,1:end-1);
    V     = Z(:,n+1:2*n);
    dvel  = V(:,2:end) - V(:,1:end-1);    % Nt x Ns
    p1    = Z(:,2*n+(1:Ns));
    p2    = Z(:,2*n+Ns+(1:Ns));
    Q     = Z(:,2*n+2*Ns+(1:Ns));
    T_o   = Z(:,end-1);

    nd = 1; if isfield(hyd,'n_parallel'), nd = max(1, hyd.n_parallel); end
    Ao = nd * (orf.n_orf * (pi*geom.d_o^2/4));
    Ap = nd * geom.Ap;

    mu_raw  = therm.mu_ref * exp(therm.b_mu*(T_o - therm.T_ref_C));
    if cfg.use_thermal && cfg.on.mu_floor
        mu = max(num.mu_min_phys, mu_raw);
    else
        mu = cfg.use_thermal.*mu_raw + (~cfg.use_thermal).*therm.mu_ref;
    end
    if cfg.use_thermal
        rho = max(100, therm.rho_ref ./ (1 + therm.alpha_rho*(T_o - therm.T_ref_C)));
    else
        rho = therm.rho_ref;
    end

    if cfg.use_thermal, p_vap=p_vap_Antoine(T_o,therm,orf); else, p_vap=p_vap_Antoine(therm.T_ref_C,therm,orf); end
    if cfg.on.cavitation, p2_eff = max(p2, orf.cav_sf * p_vap); else, p2_eff = p2; end

    if cfg.use_orifice
        if cfg.on.CdRe
            Re = rho .* abs(Q) .* geom.d_o ./ max(Ao .* mu, 1e-9);
            Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
            Cd = max(min(Cd,1.2),0.2);
        else
            Cd = orf.CdInf;
        end
        if cfg.on.Rkv, RQ = rho ./ max(2*(Cd.*Ao).^2,1e-12); else, RQ = 0*Q; end

        Qcap = getfield_default(num,'Qcap_big', ...
            0.4*( max(max(orf.CdInf,orf.Cd0)*Ao,1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ));
        if cfg.on.Qsat, Q_h = Qcap * tanh(Q./max(Qcap,1e-9)); else, Q_h = Q; end

        dP_kv = RQ .* Q_h .* abs(Q_h);
        if cfg.on.Rlam
            R_lam_t = (128 * mu .* geom.Lori ./ (pi * geom.d_o^4)) / max(1, nd*orf.n_orf);
            dP_lam = R_lam_t .* Q;
        else
            dP_lam = 0*Q;
        end
        dP_h = dP_lam + dP_kv;
    else
        dP_h = 0*Q;
    end

    dP_raw = p1 - p2_eff;
    if cfg.on.dP_cap && isfinite(num.dP_cap)
        dP_eff = num.dP_cap * tanh(dP_raw ./ max(num.dP_cap,1));
    else
        dP_eff = dP_raw;
    end
    epsm = max(1e3, double(num.softmin_eps));
    dP_orf_t = 0.5*( dP_eff + dP_h - sqrt((dP_eff - dP_h).^2 + epsm^2) );

    w_pf_vec = pf_weight(t,cfg) * cfg.PF.gain;

    if cfg.on.pressure_force
        dP_pf = dP_orf_t;                          % Nt x Ns
        if (isfield(cfg,'PF') && isfield(cfg.PF,'resistive_only') && cfg.PF.resistive_only) || (isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only)
            s = tanh(20*dvel);   % 20 ~ 1/(0.05 m/s) benzeri yumuÅŸatma
            dP_pf = s .* max(0, s .* dP_pf);
        end
        F_story = k_sd*drift + (w_pf_vec .* dP_pf) * Ap;
    else
        F_story = k_sd*drift;
    end

    F = zeros(n,Nt);
    F(1,:) = -F_story(:,1).';
    if n>2, F(2:n-1,:) = (F_story(:,1:end-1) - F_story(:,2:end)).'; end
    F(n,:) =  F_story(:,end).';
end
% ---- Buhar basÄ±ncÄ± (Antoine) -----------------------------------------
function p_v = p_vap_Antoine(T_C, therm, ~)
    if isfield(therm,'antoine_A') && isfield(therm,'antoine_B') && isfield(therm,'antoine_C')
        A=therm.antoine_A; B=therm.antoine_B; C=therm.antoine_C;
    else, A=5.0; B=1700; C=-80; end
    T_C=double(T_C); p_v = 10.^(A - B./(C + T_C));
p_v = min(max(p_v, 5), 5e2);     % 5â€“500 Pa

end

