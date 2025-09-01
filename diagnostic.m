%% Diagnostic metrics for orifice damper model
% This script runs the damper model with orifice and thermal effects
% enabled and computes diagnostic metrics per story. Results are written
% to out/diagnostic.csv.

clear; clc;

% --- Load earthquake acceleration input ---
S  = load('acc_matrix.mat','acc_matrix7');
t  = S.acc_matrix7(:,1);
ag = S.acc_matrix7(:,2);
[t,iu] = unique(t,'stable'); ag = ag(iu);
dt = median(diff(t));
t  = (t(1):dt:t(end)).';
ag = interp1(S.acc_matrix7(:,1), S.acc_matrix7(:,2), t, 'linear');

% --- Load structural and damper parameters ---
parametreler;

use_orifice = true;
use_thermal = true;

[x,a,diag] = mck_with_damper(t,ag,M,C0,K, k_sd,c_lam0, use_orifice, orf, rho, Ap, A_o, Qcap_big, mu_ref, ...
    use_thermal, thermal, T0_C, T_ref_C, b_mu, c_lam_min, c_lam_cap, Lgap, ...
    cp_oil, cp_steel, steel_to_oil_mass_ratio, toggle_gain, story_mask, ...
    n_dampers_per_story, resFactor);

% --- Per-story statistics ---
nStories = size(diag.drift,2);
q = @(x,p) quantile(x,p);

rows = (1:nStories).';
T = table(rows, ...
    q(abs(diag.drift),0.95).', ...
    q(abs(diag.story_force),0.95).', ...
    q(diag.Q,0.50).', q(diag.Q,0.95).', ...
    q(diag.dP_orf,0.50).', q(diag.dP_orf,0.95).', ...
    100*mean(diag.dP_orf<0).', ...
    q(abs(diag.PF),0.95).', ...
    repmat(diag.T_oil(end),nStories,1), ...
    repmat(diag.T_steel(end),nStories,1), ...
    repmat(diag.mu(end),nStories,1), ...
    repmat(diag.energy(end),nStories,1), ...
    'VariableNames',{ 'story','drift_p95','story_force_p95','Q_q50','Q_q95', ...
    'dP_orf_q50','dP_orf_q95','cav_pct','PF_p95','T_oil_end','T_steel_end','mu_end','energy_tot'});

if ~exist('out','dir'), mkdir('out'); end
writetable(T,'out/diagnostic.csv');

%% ---------------------------------------------------------------
%% Damper model with diagnostics
function [x,a,diag] = mck_with_damper(t,ag,M,C,K, k_sd,c_lam0, use_orf,orf,rho,Ap,Ao,Qcap, mu_ref, ...
    use_thermal, thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, toggle_gain, story_mask, ...
    n_dampers_per_story, resFactor)

    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0 = zeros(2*n,1);
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    % Story vectors
    nStories = n-1;
    Rvec = toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
    mask = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi= (mask .* ndps).';
    Rvec = Rvec.';
    Nvec = 1:nStories; Mvec = 2:n;

    % Initial temperature and viscosity
    Tser = T0_C*ones(numel(t),1);
    mu_abs = mu_ref;
    c_lam = c_lam0;

    % Solve ODE
    odef = @(tt,z) [ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - dev_force(z(1:n),z(n+1:end),c_lam,mu_abs) - M*r*agf(tt) ) ];
    sol  = ode15s(odef,[t(1) t(end)],z0,opts);
    z    = deval(sol,t).';
    x    = z(:,1:n); v = z(:,n+1:end);

    drift = x(:,Mvec) - x(:,Nvec);
    dvel  = v(:,Mvec) - v(:,Nvec);
    drift_p = drift .* Rvec; dvel_p = dvel .* Rvec;
    F_lin_p = k_sd*drift_p + c_lam*dvel_p;

    if use_orf
        qmag = Qcap * tanh( (Ap/Qcap) * sqrt(dvel_p.^2 + orf.veps^2) );
        Re   = (rho .* qmag ./ max(Ao*mu_abs,1e-9)) .* max(orf.d_o,1e-9);
        Cd   = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
        dP_calc = 0.5*rho .* ( qmag ./ max(Cd.*Ao,1e-9) ).^2;
        p_up   = orf.p_amb + abs(F_lin_p)./max(Ap,1e-12);
        dP_cav = max( (p_up - orf.p_cav_eff).*orf.cav_sf, 0 );
        dP_orf = softmin(dP_calc,dP_cav,1e5);
        sgn = dvel_p ./ sqrt(dvel_p.^2 + orf.veps^2);
        F_orf_p = dP_orf .* Ap .* sgn;
        F_p = F_lin_p + F_orf_p;
        Q = Ap * sqrt(dvel_p.^2 + orf.veps^2);
        P_orf_per = dP_orf .* Q;
    else
        F_p = F_lin_p; Q = 0*dvel_p; dP_orf = 0*dvel_p; P_orf_per = 0*dvel_p;
    end

    F_story = F_p .* (Rvec .* multi);
    P_visc_per = c_lam * (dvel_p.^2);
    P_sum = sum( (P_visc_per + P_orf_per) .* multi, 2 );
    energy = cumtrapz(t,P_sum);

    % Thermal response
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

    % Node forces for acceleration
    F = zeros(numel(t),n);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    diag = struct('drift',drift,'dvel',dvel,'story_force',F_story,'Q',Q, ...
        'dP_orf',dP_orf,'PF',F_p,'T_oil',Tser,'T_steel',Tser,'mu',mu, ...
        'energy',energy,'P_sum',P_sum,'c_lam',c_lam);

    function Fd = dev_force(x_,v_,c_lam_loc,mu_abs_loc)
        drift_ = x_(Mvec) - x_(Nvec);
        dvel_  = v_(Mvec) - v_(Nvec);
        drift_p_ = drift_ .* Rvec;
        dvel_p_  = dvel_  .* Rvec;
        F_lin_p_ = k_sd*drift_p_ + c_lam_loc*dvel_p_;
        if ~use_orf
            F_p_ = F_lin_p_;
        else
            qmag_ = Qcap * tanh( (Ap/Qcap) * sqrt(dvel_p_.^2 + orf.veps^2) );
            Re_   = (rho .* qmag_ ./ max(Ao*mu_abs_loc,1e-9)) .* max(orf.d_o,1e-9);
            Cd_   = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re_./orf.Rec).^orf.p_exp);
            dP_calc_ = 0.5*rho .* ( qmag_ ./ max(Cd_.*Ao,1e-9) ).^2;
            p_up_   = orf.p_amb + abs(F_lin_p_)./max(Ap,1e-12);
            dP_cav_ = max( (p_up_ - orf.p_cav_eff).*orf.cav_sf, 0 );
            dP_orf_ = softmin(dP_calc_,dP_cav_,1e5);
            sgn_ = dvel_p_ ./ sqrt(dvel_p_.^2 + orf.veps^2);
            F_orf_p_ = dP_orf_ .* Ap .* sgn_;
            F_p_ = F_lin_p_ + F_orf_p_;
        end
        F_story_ = F_p_ .* (Rvec .* multi);
        Fd = zeros(n,1);
        Fd(Nvec) = Fd(Nvec) - F_story_.';
        Fd(Mvec) = Fd(Mvec) + F_story_.';
    end
end

function y = softmin(a,b,epsm)
    y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
end
