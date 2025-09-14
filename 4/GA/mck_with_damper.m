function [x,a,diag] = mck_with_damper(t,ag,M,C,K, k_sd,c_lam0,Lori, use_orf,orf,rho,Ap,Ao,Qcap, mu_ref, ...
    use_thermal, thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, story_mask, ...
    n_dampers_per_story, resFactor, cfg)
%% Girdi Parametreleri
    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0 = zeros(2*n,1);
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    % Kat vektörleri
    nStories = n-1;
    mask = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi= (mask .* ndps).';
    Nvec = 1:nStories; Mvec = 2:n;

    % Başlangıç sıcaklığı ve viskozitesi
    Tser = T0_C*ones(numel(t),1);
    mu_abs = mu_ref;
    mu_min_phys = NaN;
    if isfield(cfg,'on') && isfield(cfg.on,'mu_floor') && cfg.on.mu_floor
        if isfield(cfg,'num') && isfield(cfg.num,'mu_min_phys') && isfinite(cfg.num.mu_min_phys)
            mu_min_phys = cfg.num.mu_min_phys;
        end
        mu_abs = max(mu_abs, mu_min_phys);
    end
    scale_mu = mu_abs / max(mu_ref, eps);
    c_lam = c_lam0 * scale_mu;
    R_lam_ref = c_lam0 / max(Ap^2, 1e-12);
    R_lam = R_lam_ref * scale_mu;

%% ODE Çözümü
    odef = @(tt,z) [ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - dev_force(tt,z(1:n),z(n+1:end),c_lam,mu_abs) - M*r*agf(tt) ) ];
    sol  = ode15s(odef,[t(1) t(end)],z0,opts);
    z    = deval(sol,t).';
    x    = z(:,1:n); v = z(:,n+1:end);

    drift = x(:,Mvec) - x(:,Nvec);
    dvel  = v(:,Mvec) - v(:,Nvec);
    % Faz 3: Lineer parçada sadece yay
    F_lin = k_sd*drift;

    if isfield(cfg,'compat_simple') && cfg.compat_simple
        Q = Ap * dvel;
        dP_lam = R_lam * Q;
        dP_kv  = 0*dvel;
        P_orf_per = 0*dvel;
        dP_kv_loc = 0*dvel;
        dP_cav_loc = 0*dvel;
    else
        if use_orf
            orf_loc = orf;
            params = struct('Ap_eff',Ap,'orf',orf_loc,'rho',rho,...
                            'Ao',Ao,'mu',mu_abs,'Lori',Lori);
            % Quadratic drop computed from orifice model
            [dP_lam, dP_kv_loc, Q, ~] = calc_orifice_force(dvel, params);
            % Cavitation-limited pressure drop
            p_up_loc  = orf.p_amb + abs(F_lin)./max(Ap,1e-12);
            dP_cav_loc= max( (p_up_loc - orf.p_cav_eff).*orf.cav_sf, 0 );
            dP_kv = Utils.softmin(dP_kv_loc, dP_cav_loc, eps);
            P_orf_per = dP_kv .* Q;
        else
            Q = Ap * dvel;
            dP_lam = R_lam * Q;
            dP_kv = 0*dvel;
            P_orf_per = 0*dvel;
            dP_kv_loc = 0*dvel;
            dP_cav_loc = 0*dvel;
        end
    end

    dP_resist = dP_lam + dP_kv;
    dp_pf = dP_resist;
    if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
        s = tanh(cfg.PF.resistive_slope*dvel);
        dp_pf = s .* max(0, s .* dp_pf);
    end
    w_pf_vec = Utils.pf_weight(t, cfg) * cfg.PF.gain;
    F_p = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;

    % Geometri ölçeklendirmesi R sadece montajda uygulanır
    F_story = F_p;
    P_visc_per = dP_lam .* Q;
    P_sum = sum( (P_visc_per + P_orf_per) .* multi, 2 );
    energy = cumtrapz(t,P_sum);
    P_orf_tot = sum(P_orf_per .* multi, 2);
    % Yapısal güç kat toplam kuvvetini kullanır; ekstra çarpan kullanılmaz
    P_struct_tot = sum(F_story .* dvel, 2);
    E_orifice = trapz(t, P_orf_tot);
    E_struct = trapz(t, P_struct_tot);
    P_mech = mean(P_struct_tot);

%% Termal Hesap (Phase 7: iki-düğüm diagnostik; dinamiğe geri besleme yok)
    nDtot = sum(multi);
    V_oil_per = resFactor*(Ap*(2*Lgap));
    m_oil_tot = nDtot*(rho*V_oil_per);
    m_steel_tot = steel_to_oil_mass_ratio*m_oil_tot;
    C_oil   = max(m_oil_tot*cp_oil,   eps);
    C_steel = max(m_steel_tot*cp_steel, eps);
    T_o = Tser; T_s = T0_C*ones(numel(t),1);
    hA_os   = Utils.getfield_default(thermal, 'hA_os',    thermal.hA_W_perK);
    hA_o_env= Utils.getfield_default(thermal, 'hA_o_env', thermal.hA_W_perK);
    hA_s_env= Utils.getfield_default(thermal, 'hA_s_env', thermal.hA_W_perK);
    dtv = diff(t);
    for k=1:numel(t)-1
        Pk = 0.5*(P_sum(k)+P_sum(k+1));
        dT_o = ( Pk - hA_os*(T_o(k)-T_s(k)) - hA_o_env*(T_o(k)-thermal.T_env_C) ) / C_oil;
        dT_s = ( + hA_os*(T_o(k)-T_s(k)) - hA_s_env*(T_s(k)-thermal.T_env_C) ) / C_steel;
        T_o(k+1) = T_o(k) + dtv(k)*dT_o;
        T_s(k+1) = T_s(k) + dtv(k)*dT_s;
        T_o(k+1) = min(max(T_o(k+1), T0_C), T0_C + thermal.dT_max);
        T_s(k+1) = min(max(T_s(k+1), T0_C), T0_C + thermal.dT_max);
    end
    mu = mu_ref*exp(b_mu*(T_o - T_ref_C));
    if isfinite(mu_min_phys)
        mu = max(mu, mu_min_phys);
    end

%% Çıktı Hesabı
    % İvme için düğüm kuvvetleri
    F = zeros(numel(t),n);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    diag = struct('drift',drift,'dvel',dvel,'story_force',F_story,'Q',Q, ...
        'dP_resist',dP_resist,'dP_lam',dP_lam,'dP_orf',dP_resist,'PF',F_p,'T_oil',T_o,'T_steel',T_s,'mu',mu, ...
        'energy',energy,'P_sum',P_sum,'c_lam',c_lam,'Lori',Lori, ...
        'P_orf_per',P_orf_per,'dP_kv',dP_kv,'dP_cav',dP_cav_loc, ...
        'E_orifice',E_orifice,'E_struct',E_struct);

%% İç Fonksiyonlar
    function Fd = dev_force(tt,x_,v_,c_lam_loc,mu_abs_loc)
        drift_ = x_(Mvec) - x_(Nvec);
        dvel_  = v_(Mvec) - v_(Nvec);
        F_lin_ = k_sd*drift_;
        if isfield(cfg,'compat_simple') && cfg.compat_simple
            Q_ = Ap * dvel_;
            R_lam_ = c_lam_loc / max(Ap^2,1e-12);
            dP_resist_ = R_lam_ * Q_;
            dp_pf_ = dP_resist_;
        else
            if ~use_orf
                Q_ = Ap * dvel_;
                R_lam_ = c_lam_loc / max(Ap^2,1e-12);
                dP_resist_ = R_lam_ * Q_;
            else
                params = struct('Ap_eff',Ap,'Ao',Ao,'orf',orf,'rho',rho,...
                                'mu',mu_abs_loc,'Lori',Lori);
                [dP_lam_, dP_kv_, Q_, ~] = calc_orifice_force(dvel_, params);
                p_up_loc_  = orf.p_amb + abs(F_lin_)./max(Ap,1e-12);
                dP_cav_loc_= max( (p_up_loc_ - orf.p_cav_eff).*orf.cav_sf, 0 );
                dP_kv_ = Utils.softmin(dP_kv_, dP_cav_loc_, eps);
                dP_resist_ = dP_lam_ + dP_kv_;
            end
            dp_pf_ = dP_resist_;
            if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
                s = tanh(cfg.PF.resistive_slope*dvel_);
                dp_pf_ = s .* max(0, s .* dp_pf_);
            end
        end
        w_pf = Utils.pf_weight(tt,cfg) * cfg.PF.gain;
        F_p_ = k_sd*drift_ + (w_pf .* dp_pf_) * Ap;
        F_story_ = F_p_;
        Fd = zeros(n,1);
        Fd(Nvec) = Fd(Nvec) - F_story_;
        Fd(Mvec) = Fd(Mvec) + F_story_;
    end
    function [dP_lam, dP_kv, Q, P_orf_per] = calc_orifice_force(dvel, params)
        % Laminar and quadratic pressure drops based on Q = Ap_eff*dvel.

        Q = params.Ap_eff * dvel;

        d_o = max(params.orf.d_o, 1e-12);
        n_orf = max(params.Ao / (pi*d_o^2/4), 1);
        R_lam = (128 * params.mu * params.Lori) / (pi * d_o^4 * n_orf);
        dP_lam = R_lam * Q;

        Re   = (params.rho .* abs(Q) .* d_o) ./ max(params.Ao*params.mu,1e-12);
        Cd0   = params.orf.Cd0;
        CdInf = params.orf.CdInf;
        p_exp = params.orf.p_exp;
        Rec   = params.orf.Rec;
        Cd    = CdInf - (CdInf - Cd0) ./ (1 + (Re./max(Rec,1)).^p_exp);
        Cd    = max(min(Cd, 1.2), 0.2);

        dP_kv  = params.rho ./ (2*(Cd.*params.Ao).^2) .* Q .* abs(Q);
        P_orf_per = dP_kv .* Q;
    end
end


