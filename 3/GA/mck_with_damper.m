function [x,a,diag] = mck_with_damper(t,ag,M,C,K, k_sd,c_lam0,Lori, use_orf,orf,rho,Ap,Ao,Qcap, mu_ref, ...
    use_thermal, thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, toggle_gain, story_mask, ...
    n_dampers_per_story, resFactor, cfg, F_story_target)
%% Girdi Parametreleri
    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0 = zeros(2*n,1);
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    % Kat vektörleri
    if nargin < 32 || isempty(F_story_target), F_story_target = []; end
    nStories = n-1;
    Rvec = toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
    mask = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi= (mask .* ndps).';
    Rvec = Rvec.';
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
    F_lin = k_sd*drift + c_lam*dvel;

    if use_orf
        params = struct('Ap',Ap,'Qcap',Qcap,'orf',orf,'rho',rho,...
                        'Ao',Ao,'mu',mu_abs,'F_lin',F_lin);
        [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force(dvel, params);
        F_p = F_lin + F_orf;
    else
        F_p = F_lin; Q = 0*dvel; dP_orf = 0*dvel; P_orf_per = 0*dvel;
    end

    dp_pf = (c_lam*dvel + (F_p - k_sd*drift)) ./ Ap;
    if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
        s = tanh(20*dvel);
        dp_pf = s .* max(0, s .* dp_pf);
    end
    w_pf_vec = Utils.pf_weight(t, cfg) * cfg.PF.gain;
    F_p = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;

    % Geometri ölçeklendirmesi R sadece montajda uygulanır
    F_story = F_p .* (Rvec .* multi);
    F_story_err = story_force_error(F_story, F_story_target);
    P_visc_per = c_lam * (dvel.^2);
    P_sum = sum( (P_visc_per + P_orf_per) .* multi, 2 );
    energy = cumtrapz(t,P_sum);
    P_orf_tot = sum(P_orf_per .* multi, 2);
    % Yapısal güç kat toplam kuvvetini kullanır; ekstra çarpan kullanılmaz
    P_struct_tot = sum(F_story .* dvel, 2);
    E_orifice = trapz(t, P_orf_tot);
    E_struct = trapz(t, P_struct_tot);
    P_mech = mean(P_struct_tot);

%% Termal Hesap
    nDtot = sum(multi);
    V_oil_per = resFactor*(Ap*(2*Lgap));
    m_oil_tot = nDtot*(rho*V_oil_per);
    m_steel_tot = steel_to_oil_mass_ratio*m_oil_tot;
    C_th = max(m_oil_tot*cp_oil + m_steel_tot*cp_steel, eps);
    dtv = diff(t);
    for k=1:numel(t)-1
        Pk = 0.5*(P_sum(k)+P_sum(k+1));
        Tser(k+1) = Tser(k) + dtv(k)*( Pk/C_th - (thermal.hA_W_perK/C_th)*(Tser(k)-thermal.T_env_C) );
        Tser(k+1) = min(max(Tser(k+1), T0_C), T0_C + thermal.dT_max);
    end
    mu = mu_ref*exp(b_mu*(Tser - T_ref_C));

%% Çıktı Hesabı
    % İvme için düğüm kuvvetleri
    F = zeros(numel(t),n);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    diag = struct('drift',drift,'dvel',dvel,'story_force',F_story,'Q',Q, ...
        'dP_orf',dP_orf,'PF',F_p,'T_oil',Tser,'mu',mu, ...
        'energy',energy,'P_sum',P_sum,'c_lam',c_lam,'Lori',Lori, ...
        'E_orifice',E_orifice,'E_struct',E_struct,'F_story_err',F_story_err);
    if ~isempty(F_story_target)
        diag.F_story_target = F_story_target;
    end

%% İç Fonksiyonlar
    function Fd = dev_force(tt,x_,v_,c_lam_loc,mu_abs_loc)
        Rcol = Rvec.';
        multicol = multi.';
        drift_ = x_(Mvec) - x_(Nvec);
        dvel_  = v_(Mvec) - v_(Nvec);
        F_lin_ = k_sd*drift_ + c_lam_loc*dvel_;
        if ~use_orf
            F_orf_ = 0*dvel_;
        else
            params = struct('Ap',Ap,'Qcap',Qcap,'orf',orf,'rho',rho,...
                            'Ao',Ao,'mu',mu_abs_loc,'F_lin',F_lin_);
            [F_orf_, ~, ~, ~] = calc_orifice_force(dvel_, params);
        end
        dp_pf_ = (c_lam_loc*dvel_ + F_orf_) ./ Ap;
        if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
            s = tanh(20*dvel_);
            dp_pf_ = s .* max(0, s .* dp_pf_);
        end
        w_pf = Utils.pf_weight(tt,cfg) * cfg.PF.gain;
        F_p_ = k_sd*drift_ + (w_pf .* dp_pf_) * Ap;
        % R ölçeklendirmesi yalnızca montajda uygulanır (R*multi)
        F_story_ = F_p_ .* (Rcol .* multicol);
        Fd = zeros(n,1);
        Fd(Nvec) = Fd(Nvec) - F_story_;
        Fd(Mvec) = Fd(Mvec) + F_story_;
    end
    function [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force(dvel, params)
        qmag = params.Qcap * tanh( (params.Ap/params.Qcap) * sqrt(dvel.^2 + params.orf.veps^2) );
        Re   = (params.rho .* qmag ./ max(params.Ao*params.mu,1e-9)) .* max(params.orf.d_o,1e-9);
        Cd0   = params.orf.Cd0;
        CdInf = params.orf.CdInf;
        p_exp = params.orf.p_exp;
        Cd   = CdInf - (CdInf - Cd0) ./ (1 + (Re./params.orf.Rec).^p_exp);
        dP_calc = 0.5*params.rho .* ( qmag ./ max(Cd.*params.Ao,1e-9) ).^2;
        p_up   = params.orf.p_amb + abs(params.F_lin)./max(params.Ap,1e-12);
        dP_cav = max( (p_up - params.orf.p_cav_eff).*params.orf.cav_sf, 0 );
        dP_orf = Utils.softmin(dP_calc,dP_cav,1e5);
        sgn = dvel ./ sqrt(dvel.^2 + params.orf.veps^2);
        F_orf = dP_orf .* params.Ap .* sgn;
        Q = params.Ap * sqrt(dvel.^2 + params.orf.veps^2);
        P_orf_per = dP_orf .* Q;
    end
    function err = story_force_error(F_actual, F_target)
        if nargin < 2 || isempty(F_target)
            err = 0;
            return;
        end
        n = min(size(F_actual,1), size(F_target,1));
        diff = F_actual(1:n,:) - F_target(1:n,:);
        err = sqrt(mean(diff(:).^2));
    end
end
