%% ================================================================
%  10-Katlı Çerçeve — ODE-only (ham ivme)
%  Dampersiz ve damperli (k_sd + c_lam + ORIFIS) üst üste
%  NOT: Damper hikâye kuvvetleri düğümlere F = G^T*q ile dağıtılır
% ================================================================

clear; clc; close all;

%% 0) Deprem girdisi (ham ivme, m/s^2)
S = load('acc_matrix.mat','acc_matrix');
t  = S.acc_matrix(:,1);
ag = S.acc_matrix(:,2);

% tekilleştir + eş-adımlı ufak düzeltme
[t,iu] = unique(t,'stable'); ag = ag(iu);
dt0 = median(diff(t)); 
t  = (t(1):dt0:t(end)).';
ag = interp1(S.acc_matrix(:,1), S.acc_matrix(:,2), t, 'linear');

% Arias %5–%95 (sadece gösterim)
[t5,t95] = arias_win(t,ag,0.05,0.95);

%% 1) Yapı (10 kat)
n  = 10;
m  = 360e3 * ones(n,1);
k  = 6.5e8 * ones(n,1);
c  = 6.2e6 * ones(n,1);
[M,K,C0] = make_KCM(n,m,k,c);

% T1 (referans)
[~,D] = eig(K,M); w = sqrt(sort(diag(D),'ascend')); T1 = 2*pi/w(1);

%% 2) Damper parametreleri (senin verdiğin)
Dp=0.125;   Lgap=0.055;  d_o=2.7e-3;  Lori=0.10;  mu=0.9;     % [Pa·s]
Kd=1.6e9;   Ebody=2.1e11; Gsh=79e9;   d_w=12e-3;  D_m=80e-3; n_turn=8;

Ap   = pi*Dp^2/4;
k_h  = Kd*Ap^2/Lgap;
k_s  = Ebody*Ap/Lgap;
k_hyd= 1/(1/k_h + 1/k_s);
k_p  = Gsh*d_w^4/(8*n_turn*D_m^3);
k_sd = k_hyd + k_p;
c_lam= 12*mu*Lori*Ap^2/d_o^4;

% Orifis (tam model – Cd(Re) + soft-min kavitasyon)
rho        = 850;                  % yağ yoğunluğu [kg/m^3]
orf.Cd0    = 0.60;
orf.CdInf  = 0.82;
orf.Rec    = 3e3;
orf.p_exp  = 1.2;
orf.p_amb  = 1e5;                  % ortam basıncı [Pa]
orf.p_vap  = 3e3;                  % buhar basıncı [Pa]
orf.cav_sf = 0.85;                 % güvenlik kats.
orf.n_orf  = 2;                    % orifis adedi
geom.Ap    = Ap;
geom.Ao    = orf.n_orf*(pi*d_o^2/4);
orf.enable = true;   % <- ORIFIS AÇ/KAPA. false yaparsan saf (k_sd + c_lam) modele döner.


%% 3) ODE çözümleri
[x0,~]  = lin_MCK(t,ag,M,C0,K);                                       % dampersiz
[xd,~]  = mck_with_damper_orf(t,ag,M,C0,K,k_sd,c_lam,geom,orf,rho,mu);% damper+orifis

x10_0 = x0(:,10);
x10_d = xd(:,10);

%% 4) Grafikler
figure('Name','x10(t) — ham ivme (ODE-only)','Color','w');
plot(t,x10_0,'b','LineWidth',1.4); hold on;
plot(t,x10_d,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off'); plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('x_{10}(t) [m]');
title(sprintf('10. Kat yer değiştirme — ham ivme (ODE-only) | T1=%.3f s',T1));
legend('dampersiz','damperli (orifisli)','Location','best');

% Arias penceresi içinde detrend karşılaştırma
idx = (t>=t5 & t<=t95);
tt  = t(idx);
d0  = detrend(x10_0(idx),'linear');
dd  = detrend(x10_d(idx),'linear');
figure('Name','Arias penceresi — detrended x10','Color','w');
plot(tt,d0,'b','LineWidth',1.2); hold on;
plot(tt,dd,'r','LineWidth',1.0); grid on;
xlabel('t [s]'); ylabel('x_{10} detrended [m]'); title('Arias penceresi');

fprintf('x10_max  (dampersiz)          = %.4g m\n',max(abs(x10_0)));
fprintf('x10_max  (damperli + orifis) = %.4g m\n',max(abs(x10_d)));

%% ===================== Yardımcı fonksiyonlar =====================

function [M,K,C] = make_KCM(n,mv,kv,cv)
    M = diag(mv);
    K = zeros(n); C = zeros(n);
    for i=1:n
        kL=kv(i); cL=cv(i);
        if i<n, kU=kv(i+1); cU=cv(i+1); else, kU=0; cU=0; end
        K(i,i) = kL + (i<n)*kU;    C(i,i) = cL + (i<n)*cU;
        if i>1, K(i,i-1)=-kL; C(i,i-1)=-cL; end
        if i<n, K(i,i+1)=-kU; C(i,i+1)=-cU; end
    end
end

function [x,a] = lin_MCK(t,ag,M,C,K)
    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    odef = @(tt,z)[ z(n+1:end);
        M \ ( -C*z(n+1:end) - K*z(1:n) - M*r*agf(tt) ) ];
    z0 = zeros(2*n,1);
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    sol = ode15s(odef,[t(1) t(end)],z0,opts);
    z = deval(sol,t).';
    x = z(:,1:n);
    a = ( -(M\(C*z(:,n+1:end).' + K*z(:,1:n).')).' - ag.*r.' );
end

function [x,a] = mck_with_damper_orf(t,ag,M,C,K,k_sd,c_lam,geom,orf,rho,mu)
    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    odef = @(tt,z) rhs(tt,z);
    z0 = zeros(2*n,1);
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    sol = ode15s(odef,[t(1) t(end)],z0,opts);
    z = deval(sol,t).';
    x = z(:,1:n); v = z(:,n+1:end);
    Fdev = dev_force_orf(x,v,k_sd,c_lam,geom,orf,rho,mu);
    a = ( -(M\(C*v.' + K*x.' + Fdev)).' - ag.*r.' );

    function dz = rhs(tt,zz)
        x_ = zz(1:n); v_ = zz(n+1:end);
        Fd = dev_force_orf(x_,v_,k_sd,c_lam,geom,orf,rho,mu);
        dv = M \ ( -C*v_ - K*x_ - Fd - M*r*agf(tt) );
        dz = [v_; dv];
    end
end

function F = dev_force_orf(x,v,k_sd,c_lam,geom,orf,rho,mu)
    % x,v: (n x 1) veya (Nt x n). Çıkış F: düğüm kuvvetleri (n x 1) veya (n x Nt)
    if isvector(x), x=x(:).'; v=v(:).'; end
    Nt = size(x,1); n = size(x,2);

    % Kat farkları
    drift = x(:,2:end) - x(:,1:end-1);
    dvel  = v(:,2:end) - v(:,1:end-1);

    % Lineer parça
    F_lin = k_sd.*drift + c_lam.*dvel;

    % --- Orifis (Cd(Re) + soft-min kavitasyon) ---
    Ap = geom.Ap; Ao = geom.Ao;

    % numerik yumuşatma ve debi satürasyonu
    Veps = 0.20;                                  % |dvel| için yumuşatma
    Qcap = max(orf.CdInf*Ao, 1e-9) * sqrt(2*1e9/rho);
    qmag = Qcap * tanh((Ap/Qcap)*sqrt(dvel.^2 + Veps^2));   % ~Ap*|dvel| satürasyonlu

    Re    = rho .* qmag ./ max(Ao*mu, 5e-9);      % ~ρ*V*d/μ (d sabit → Rec'a gömülü)
    CdEff = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);

    dP_calc = 0.5*rho .* ( qmag ./ max(CdEff.*Ao, 5e-9) ).^2;

    % kavitasyon üst sınırı (soft-min)
    p_up   = orf.p_amb + abs(F_lin)./max(Ap,1e-9);
    dP_cav = max((p_up - orf.p_vap).*orf.cav_sf, 0);
    dP_orf = softmin(dP_calc, dP_cav, 1e5);       % 1e5 Pa yumuşatma ölçeği

    sgn    = dvel ./ sqrt(dvel.^2 + Veps^2);
    F_orf  = dP_orf .* Ap .* sgn;

    F_story = F_lin + F_orf;

    % Düğümlere dağıt (F = G^T q)  —  (önceki çalışan işaretle aynı)
    F = zeros(Nt,n);
    F(:,1)      = -F_story(:,1);
    if n>2, F(:,2:n-1) = F_story(:,1:end-1) - F_story(:,2:end); end
    F(:,n)      =  F_story(:,end);
    F = F.';  if size(F,2)==1, F = F(:,1); end
end

function y = softmin(a,b,epsm)
    % C2-sürekli yaklaşık min(a,b)
    y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
end

function [t5,t95] = arias_win(t,ag,p1,p2)
    if nargin<3, p1=0.05; p2=0.95; end
    a2 = ag(:).^2; dt = [diff(t); median(diff(t))];
    E  = cumsum(0.5*(a2+[a2(2:end);a2(end)]).*dt);
    t5  = interp1(E,t,p1*E(end),'linear','extrap');
    t95 = interp1(E,t,p2*E(end),'linear','extrap');
end
