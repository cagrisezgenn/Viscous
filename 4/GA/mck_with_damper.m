function [x,a_rel,ts] = mck_with_damper(t,ag,M,C,K, k_sd,c_lam0,Lori, orf,rho,Ap,Ao,Qcap, mu_ref, ...
    thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
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
    c_lam = c_lam0;

%% ODE Çözümü
    odef = @(tt,z) [ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - dev_force(tt,z(1:n),z(n+1:end),c_lam,mu_abs) - M*r*agf(tt) ) ];
    sol  = ode15s(odef,[t(1) t(end)],z0,opts);
    z    = deval(sol,t).';
    x    = z(:,1:n); v = z(:,n+1:end);

    drift = x(:,Mvec) - x(:,Nvec);
    dvel  = v(:,Mvec) - v(:,Nvec);
    % Faz 3: Lineer parçada sadece yay (laminer PF tarafında)
    F_lin = k_sd*drift;

    % Faz 6: Qcap ölçeği ve softmin eps opsiyonu
    Qcap_eff = Qcap;
    if isfield(cfg,'num') && isfield(cfg.num,'Qcap_scale') && isfinite(cfg.num.Qcap_scale)
        Qcap_eff = max(1e-9, Qcap * cfg.num.Qcap_scale);
    end
    orf_loc = orf;
    if isfield(cfg,'num') && isfield(cfg.num,'softmin_eps') && isfinite(cfg.num.softmin_eps)
        orf_loc.softmin_eps = cfg.num.softmin_eps;
    end
    params = struct('Ap',Ap,'Qcap',Qcap_eff,'orf',orf_loc,'rho',rho,...
                    'Ao',Ao,'mu',mu_abs,'F_lin',F_lin,'Lori',Lori);
    [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force(dvel, params);
    % Ek diagnostikler: dP_kv ve dP_cav (kv ve kavitasyon limitleri)
    qmag_loc = Qcap_eff * tanh( (Ap/Qcap_eff) * sqrt(dvel.^2 + orf.veps^2) );
    Re_loc   = (rho .* qmag_loc .* max(orf.d_o,1e-9)) ./ max(Ao*mu_abs,1e-9);
    Cd_loc0  = orf.Cd0; Cd_locInf = orf.CdInf; Rec_loc = orf.Rec; pexp_loc = orf.p_exp;
    Cd_loc = Cd_locInf - (Cd_locInf - Cd_loc0) ./ (1 + (Re_loc./max(Rec_loc,1)).^pexp_loc);
    Cd_loc = max(min(Cd_loc,1.2),0.2);
    dP_kv_loc = 0.5*rho .* ( qmag_loc ./ max(Cd_loc.*Ao,1e-12) ).^2;
    p_up_loc  = orf.p_amb + abs(F_lin)./max(Ap,1e-12);
    dP_cav_loc= max( (p_up_loc - orf.p_cav_eff).*orf.cav_sf, 0 );
    F_p = F_lin + F_orf;

    dp_pf = (c_lam*dvel + (F_p - k_sd*drift)) ./ Ap;
    if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
        s = tanh(20*dvel);
        dp_pf = s .* max(0, s .* dp_pf);
    end
    w_pf_vec = Utils.pf_weight(t, cfg) * cfg.PF.gain;
    F_p = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;

    % Geometri ölçeklendirmesi R sadece montajda uygulanır
    F_story = F_p;
    P_visc_per = c_lam * (dvel.^2);
    P_sum = sum( (P_visc_per + P_orf_per) .* multi, 2 );
    P_orf_tot = sum(P_orf_per .* multi, 2);
    % Yapısal güç kat toplam kuvvetini kullanır; ekstra çarpan kullanılmaz
    P_struct_tot = sum(F_story .* dvel, 2);
    E_orf = cumtrapz(t, P_orf_tot);
    E_struct = cumtrapz(t, P_struct_tot);

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

%% Çıktı Hesabı
    % İvme için düğüm kuvvetleri
    F = zeros(numel(t),n);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a_rel = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    ts = struct();
    ts.dvel = dvel;
    ts.story_force = F_story;
    ts.Q = Q;
    ts.dP_orf = dP_orf;
    ts.PF = F_p;
    ts.cav_mask = dP_orf < 0;
    ts.P_sum = P_sum;
    ts.E_orf = E_orf;
    ts.E_struct = E_struct;
    ts.T_oil = T_o;
    ts.mu = mu;
    ts.c_lam = c_lam;

%% İç Fonksiyonlar
    function Fd = dev_force(tt,x_,v_,c_lam_loc,mu_abs_loc)
        drift_ = x_(Mvec) - x_(Nvec);
        dvel_  = v_(Mvec) - v_(Nvec);
        % Sütun yönelimli etkin parametreler
        % Faz 3: Lineer parçada sadece yay
        F_lin_ = k_sd*drift_;
        params = struct('Ap',Ap,'Qcap',Qcap,'orf',orf,'rho',rho,...
                        'Ao',Ao,'mu',mu_abs_loc,'F_lin',F_lin_,'Lori',Lori);
        [F_orf_, ~, ~, ~] = calc_orifice_force(dvel_, params);
        dp_pf_ = (c_lam_loc*dvel_ + F_orf_) ./ Ap;
        if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
            s = tanh(20*dvel_);
            dp_pf_ = s .* max(0, s .* dp_pf_);
        end
        w_pf = Utils.pf_weight(tt,cfg) * cfg.PF.gain;
        F_p_ = k_sd*drift_ + (w_pf .* dp_pf_) * Ap;
        % R ölçeklendirmesi yalnızca montajda uygulanır (R*multi)
        F_story_ = F_p_;
        Fd = zeros(n,1);
        Fd(Nvec) = Fd(Nvec) - F_story_;
        Fd(Mvec) = Fd(Mvec) + F_story_;
    end
    function [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force(dvel, params)
        % Phase 6 (no p-states): smoother Cd(Re) and kv-only orifice drop.
        % Laminar viscous loss is accounted in PF via c_lam*dvel; avoid double counting here.

        % Saturated volumetric flow magnitude (stability)
        qmag = params.Qcap * tanh( (params.Ap/params.Qcap) * sqrt(dvel.^2 + params.orf.veps^2) );

        % Reynolds and discharge coefficient (clamped)
        Re   = (params.rho .* qmag .* max(params.orf.d_o,1e-9)) ./ max(params.Ao*params.mu,1e-9);
        Cd0   = params.orf.Cd0;
        CdInf = params.orf.CdInf;
        p_exp = params.orf.p_exp;
        Rec   = params.orf.Rec;
        Cd    = CdInf - (CdInf - Cd0) ./ (1 + (Re./max(Rec,1)).^p_exp);
        Cd    = max(min(Cd, 1.2), 0.2);

        % kv-only drop
        dP_kv  = 0.5*params.rho .* ( qmag ./ max(Cd.*params.Ao,1e-12) ).^2;

        % Cavitation soft-limit via softmin
        p_up   = params.orf.p_amb + abs(params.F_lin)./max(params.Ap,1e-12);
        dP_cav = max( (p_up - params.orf.p_cav_eff).*params.orf.cav_sf, 0 );
        epsm = 1e5;
        if isfield(params,'orf') && isfield(params.orf,'softmin_eps') && isfinite(params.orf.softmin_eps)
            epsm = params.orf.softmin_eps;
        end
        dP_orf = Utils.softmin(dP_kv, dP_cav, epsm);

        % Force sign from velocity (no p-states)
        sgn = dvel ./ sqrt(dvel.^2 + params.orf.veps^2);
        F_orf = dP_orf .* params.Ap .* sgn;

        % Diagnostics (positive)
        Q = qmag;
        P_orf_per = dP_kv .* qmag;   % avoid counting laminar twice
    end
end


