function [x,a,diag] = mck_with_damper(t,ag,M,C,K, k_sd,c_lam0,Lori, use_orf,orf,rho,Ap,Ao,Qcap, mu_ref, ...
    use_thermal, thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, toggle_gain, story_mask, ...
    n_dampers_per_story, resFactor, cfg)
%% Girdi Parametreleri
    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0 = zeros(2*n+1,1); z0(end) = T0_C; % son eleman sıcaklık
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    % Kat vektörleri
    nStories = n-1;
    Rvec = toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
    mask = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi= (mask .* ndps);      % kolon vektör
    Rcol = Rvec(:);             % R ölçeği kolon
    Nvec = 1:nStories; Mvec = 2:n;

    % Termal kapasite
    nDtot = sum(multi);
    V_oil_per = resFactor*(Ap*(2*Lgap));
    m_oil_tot = nDtot*(rho*V_oil_per);
    m_steel_tot = steel_to_oil_mass_ratio*m_oil_tot;
    C_th = max(m_oil_tot*cp_oil + m_steel_tot*cp_steel, eps);

%% ODE Çözümü
    odef = @(tt,z) local_rhs(tt,z);
    sol  = ode15s(odef,[t(1) t(end)],z0,opts);
    Z    = deval(sol,t).';
    x    = Z(:,1:n); v = Z(:,n+1:2*n); Tser = Z(:,end);

    drift = x(:,Mvec) - x(:,Nvec);
    dvel  = v(:,Mvec) - v(:,Nvec);
    mu = mu_ref*exp(b_mu*(Tser - T_ref_C));
    if use_thermal
        c_lam = min(max(c_lam0*(mu/mu_ref), c_lam_min), c_lam_cap);
    else
        c_lam = c_lam0*ones(size(mu));
    end

    F_lin = k_sd*drift + c_lam.*dvel;

    if use_orf
        qmag = Qcap * tanh( (Ap/Qcap) * sqrt(dvel.^2 + orf.veps^2) );
        Re   = (rho .* qmag ./ max(Ao*mu,1e-9)) .* max(orf.d_o,1e-9);
        Cd   = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
        dP_calc = 0.5*rho .* ( qmag ./ max(Cd.*Ao,1e-9) ).^2;
        p_up   = orf.p_amb + abs(F_lin)./max(Ap,1e-12);
        dP_cav = max( (p_up - orf.p_cav_eff).*orf.cav_sf, 0 );
        dP_orf = Utils.softmin(dP_calc,dP_cav,1e5);
        sgn = dvel ./ sqrt(dvel.^2 + orf.veps^2);
        F_orf = dP_orf .* Ap .* sgn;
        Q = Ap * sqrt(dvel.^2 + orf.veps^2);
        P_orf_per = dP_orf .* Q;
    else
        F_orf = 0*dvel; Q = 0*dvel; dP_orf = 0*dvel; P_orf_per = 0*dvel;
    end

    dp_pf = (c_lam.*dvel + F_orf) ./ Ap;
    if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
        s = tanh(20*dvel);
        dp_pf = s .* max(0, s .* dp_pf);
    end
    w_pf_vec = Utils.pf_weight(t, cfg) * cfg.PF.gain;
    F_p = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;

    F_story = F_p .* (Rvec.' .* multi.');
    P_visc_per = c_lam .* (dvel.^2);
    P_sum = sum( (P_visc_per + P_orf_per) .* multi.', 2 );
    energy = cumtrapz(t,P_sum);
    P_orf_tot = sum(P_orf_per .* multi.', 2);
    P_struct_tot = sum(F_story .* dvel, 2);
    E_orifice = trapz(t, P_orf_tot);
    E_struct = trapz(t, P_struct_tot);

%% Çıktı Hesabı
    F = zeros(numel(t),n);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    diag = struct('drift',drift,'dvel',dvel,'story_force',F_story,'Q',Q, ...
        'dP_orf',dP_orf,'PF',F_p,'T_oil',Tser,'mu',mu, ...
        'energy',energy,'P_sum',P_sum,'c_lam',c_lam,'Lori',Lori, ...
        'E_orifice',E_orifice,'E_struct',E_struct);

%% İç Fonksiyon
    function dz = local_rhs(tt,z)
        x_ = z(1:n); v_ = z(n+1:2*n); T_ = z(end);
        mu_loc = mu_ref * exp(b_mu*(T_ - T_ref_C));
        if use_thermal
            c_lam_loc = min(max(c_lam0*(mu_loc/mu_ref), c_lam_min), c_lam_cap);
        else
            c_lam_loc = c_lam0;
        end
        drift_ = x_(Mvec) - x_(Nvec);
        dvel_  = v_(Mvec) - v_(Nvec);
        F_lin_ = k_sd*drift_ + c_lam_loc*dvel_;
        if use_orf
            qmag_ = Qcap * tanh( (Ap/Qcap) * sqrt(dvel_.^2 + orf.veps^2) );
            Re_   = (rho .* qmag_ ./ max(Ao*mu_loc,1e-9)) .* max(orf.d_o,1e-9);
            Cd_   = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re_./orf.Rec).^orf.p_exp);
            dP_calc_ = 0.5*rho .* ( qmag_ ./ max(Cd_.*Ao,1e-9) ).^2;
            p_up_   = orf.p_amb + abs(F_lin_)./max(Ap,1e-12);
            dP_cav_ = max( (p_up_ - orf.p_cav_eff).*orf.cav_sf, 0 );
            dP_orf_ = Utils.softmin(dP_calc_,dP_cav_,1e5);
            sgn_ = dvel_ ./ sqrt(dvel_.^2 + orf.veps^2);
            F_orf_ = dP_orf_ .* Ap .* sgn_;
            Q_loc = Ap * sqrt(dvel_.^2 + orf.veps^2);
            P_orf_per_ = dP_orf_ .* Q_loc;
        else
            F_orf_ = 0*dvel_;
            P_orf_per_ = 0*dvel_;
        end
        dp_pf_ = (c_lam_loc*dvel_ + F_orf_) ./ Ap;
        if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
            s = tanh(20*dvel_);
            dp_pf_ = s .* max(0, s .* dp_pf_);
        end
        w_pf = Utils.pf_weight(tt,cfg) * cfg.PF.gain;
        F_p_ = k_sd*drift_ + (w_pf .* dp_pf_) * Ap;
        F_story_ = F_p_ .* (Rcol .* multi);
        P_visc_per_ = c_lam_loc * (dvel_.^2);
        P_sum_ = sum( (P_visc_per_ + P_orf_per_) .* multi );
        Fd = zeros(n,1);
        Fd(Nvec) = Fd(Nvec) - F_story_;
        Fd(Mvec) = Fd(Mvec) + F_story_;
        a_ = M \ ( -C*v_ - K*x_ - Fd - M*r*agf(tt) );
        if use_thermal
            dT_ = P_sum_/C_th - (thermal.hA_W_perK/C_th)*(T_ - thermal.T_env_C);
        else
            dT_ = 0;
        end
        dz = [v_; a_; dT_];
    end
end

