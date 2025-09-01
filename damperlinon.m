%% ================================================================
%  10-Katlı Çerçeve — ODE-only (ham ivme, m/s^2)
%  Dampersiz vs. Lineer Damper vs. Orifisli Damper (+ termal döngü anahtarı)
%  Notlar:
%   - use_orifice=true  → Cd(Re)+kavitasyon (soft-min) aktif
%   - use_orifice=false → lineer (k_sd + c_lam) model
%   - use_thermal=true  → ΔT relaxation + clamp + c_lam(T) alt/üst sınır
%   - Re = (rho*q/(A_o*mu))*d_o   (d_o düzeltmesi)
% ================================================================

clear; clc; close all;

%% --- Model anahtarları ---
use_orifice = tru;     % Orifis modeli aç/kapa
use_thermal = true;     % Termal döngü (ΔT ve c_lam(T)) aç/kapa

%% 0) Deprem girdisi (ham ivme, m/s^2)
S  = load('acc_matrix.mat','acc_matrix7');   % gerekirse path'i değiştirin
t  = S.acc_matrix7(:,1);
ag = S.acc_matrix7(:,2);

% tekilleştir + eş-adımlı küçük düzeltme
[t,iu] = unique(t,'stable'); ag = ag(iu);
dt  = median(diff(t));
t   = (t(1):dt:t(end)).';
ag  = interp1(S.acc_matrix7(:,1), S.acc_matrix7(:,2), t, 'linear');

% (Sadece gösterim) Arias %5–%95 penceresi
[t5,t95] = arias_win(t,ag,0.05,0.95);

%% 1) Yapı (10 kat)
n  = 10;
m  = 360e3 * ones(n,1);
k  = 6.5e8 * ones(n,1);
c  = 6.2e6 * ones(n,1);
[M,K,C0] = make_KCM(n,m,k,c);

% T1 (dampersiz referans)
[~,D] = eig(K,M); w = sqrt(sort(diag(D),'ascend')); T1 = 2*pi/w(1);

%% 2) Damper geometrisi ve malzeme (verdiğin değerler)
Dp=0.125;   Lgap=0.055;  d_o=2.7e-3;  Lori=0.10;  mu_ref=0.9;   % mu_ref: referans viskozite [Pa·s]
Kd=1.6e9;   Ebody=2.1e11; Gsh=79e9;   d_w=12e-3;  D_m=80e-3; n_turn=8;

% Türetilen sabitler (lineer eşdeğer)
Ap    = pi*Dp^2/4;
k_h   = Kd*Ap^2/Lgap;
k_s   = Ebody*Ap/Lgap;
k_hyd = 1/(1/k_h + 1/k_s);
k_p   = Gsh*d_w^4/(8*n_turn*D_m^3);
k_sd  = k_hyd + k_p;                          % (tek değer → tüm kat aralarına)
c_lam0= 12*mu_ref*Lori*Ap^2/d_o^4;            % T0'da laminer eşdeğer sönüm

%% 3) Orifis + termal parametreleri (önerilen set)
rho   = 850;           % yağ yoğunluğu [kg/m^3]
n_orf = 2;             % kat başına orifis sayısı
A_o   = n_orf * (pi*d_o^2/4);

orf.Cd0   = 0.61;
orf.CdInf = 0.80;
orf.Rec   = 3000;
orf.p_exp = 1.1;
orf.p_amb = 1.0e5;         % ortam basıncı
orf.p_cav_eff = 2.0e3;     % etkin kavitasyon eşiği (Pa)
orf.cav_sf    = 0.90;      % emniyet katsayısı
orf.d_o   = d_o;           % Re düzeltmesi için çap
orf.veps  = 0.10;          % [m/s] küçük hız yumuşatma

% Akış satürasyonu (sayısal kararlılık için)
Qcap_big = max(orf.CdInf*A_o, 1e-9) * sqrt(2*1.0e9/rho);

% --- Termal model (ΔT relaxation + clamp + c_lam(T) sınırları) ---
T0_C      = 25;       T_ref_C = 25;   b_mu = -0.013;      % μ(T)=μ_ref*exp(b_mu*(T-Tref))
thermal.hA_W_perK = 150;                                          % konvektif ısı kaybı
thermal.T_env_C   = 25;
thermal.max_iter  = 3;
thermal.tol_K     = 0.5;
thermal.relax     = 0.5;
thermal.dT_max    = 80;                                           % ΔT clamp (sim. üst sınır)
steel_to_oil_mass_ratio = 1.5;
n_dampers_per_story    = 1;
cp_oil   = 1800;   cp_steel = 500;   resFactor = 3;               % hacim/kapasite ölçek

% c_lam(T) sınırları
c_lam_cap      = 2e7;                 % üst sınır (cap)
c_lam_min_frac = 0.05;                % taban → 0.05*c_lam0
c_lam_min_abs  = 1e5;                 % mutlak taban
c_lam_min      = max(c_lam_min_abs, c_lam_min_frac*c_lam0);

%% 4) Çözümler
[x0,~]    = lin_MCK(t,ag,M,C0,K);  % dampersiz

% Lineer (termal kapalı tutarak)
[x_lin,~] = mck_with_damper( ...
    t,ag,M,C0,K, k_sd, c_lam0, false, orf, rho, Ap, A_o, Qcap_big, mu_ref, ...
    false, thermal, T0_C, T_ref_C, b_mu, c_lam_min, c_lam_cap, Lgap, ...
    cp_oil, cp_steel, steel_to_oil_mass_ratio, n_dampers_per_story, resFactor);

% Orifisli (+ termal anahtarına göre)
[x_orf,diag] = mck_with_damper( ...
    t,ag,M,C0,K, k_sd, c_lam0, use_orifice, orf, rho, Ap, A_o, Qcap_big, mu_ref, ...
    use_thermal, thermal, T0_C, T_ref_C, b_mu, c_lam_min, c_lam_cap, Lgap, ...
    cp_oil, cp_steel, steel_to_oil_mass_ratio, n_dampers_per_story, resFactor);

x10_0   = x0(:,10);
x10_lin = x_lin(:,10);
x10_orf = x_orf(:,10);

%% 5) Grafikler (üst üste)
figure('Name','10. Kat yer değiştirme — ham ivme (ODE-only)','Color','w');
plot(t,x10_0 ,'k','LineWidth',1.4); hold on;
plot(t,x10_lin,'b','LineWidth',1.1);
plot(t,x10_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('x_{10}(t) [m]');
title(sprintf('10-Kat | T1=%.3f s | Arias [%.3f, %.3f] s',T1,t5,t95));
legend('Dampersiz','Lineer damper','Orifisli damper','Location','best');

% Kısa özet
fprintf('x10_max  (dampersiz)   = %.4g m\n', max(abs(x10_0)));
fprintf('x10_max  (lineer)      = %.4g m\n', max(abs(x10_lin)));
fprintf('x10_max  (orifisli%s)  = %.4g m\n', tern(use_thermal,'+termal',''), max(abs(x10_orf)));
if use_orifice && use_thermal && isfield(diag,'dT_est')
    fprintf('Termal döngü: ΔT_est=%.2f K | c_lam(final)=%.3e N·s/m\n', diag.dT_est, diag.c_lam);
end
fprintf('Not: orifis modelini kapatmak için use_orifice=false; termali kapatmak için use_thermal=false.\n');

%% ===================== Yardımcı fonksiyonlar =====================
function [M,K,C] = make_KCM(n,mv,kv,cv)
    M = diag(mv);
    K = zeros(n); C = zeros(n);
    for i=1:n
        kL=kv(i); cL=cv(i);
        if i<n, kU=kv(i+1); cU=cv(i+1); else, kU=0; cU=0; end
        K(i,i) = kL + (i<n)*kU;   C(i,i) = cL + (i<n)*cU;
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

function [x,a,diag] = mck_with_damper( ...
    t,ag,M,C,K, k_sd,c_lam0, use_orf,orf,rho,Ap,Ao,Qcap, mu_ref, ...
    use_thermal, thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap, Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, n_dampers_per_story, resFactor)

    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0  = zeros(2*n,1);
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    % --- Termal arayıncı başlangıç ---
    T      = T0_C;              % anlık sıcaklık tahmini
    dT_est = 0;                 % toplam artış tahmini

    % Döngü: termal kapalıysa tek atım
    nIter = tern(use_thermal, thermal.max_iter, 1);

    for it = 1:nIter
        % ---- c_lam(T): μ(T)=μ_ref*exp(b_mu*(T-Tref)) ölçeği + clamp ----
        mu_scale = max(1e-6, exp(b_mu*(T - T_ref_C)));   % boyutsuz μ/μ_ref
        mu_abs   = mu_ref * mu_scale;                    % Pa·s (mutlak viskozite)
        c_raw    = c_lam0 * mu_scale;                    % c_lam(T0) * ölçek
        c_lam    = min(max(c_raw, c_lam_min), c_lam_cap);

        % ---- ODE çöz ----
        odef = @(tt,z) rhs(tt,z,c_lam,mu_abs);
        sol  = ode15s(odef,[t(1) t(end)],z0,opts);
        z    = deval(sol,t).';
        x    = z(:,1:n); v = z(:,n+1:end);

        % ---- Termal güncelle (kapalıysa atla) ----
        if ~use_thermal
            break
        end

        % Story farkları
        drift = x(:,2:end) - x(:,1:end-1);      % [Nt x (n-1)]
        dvel  = v(:,2:end) - v(:,1:end-1);

        % Viskoz güç
        P_visc = c_lam * (dvel.^2);             % [Nt x (n-1)]

        % Orifis güç (varsa)
        if use_orf
            qmag   = Qcap * tanh( (Ap/Qcap)*sqrt(dvel.^2 + orf.veps^2) );
            % Re = (rho*q/Ao/μ)*d_o  → μ_abs ile
            Re     = (rho .* qmag ./ max(Ao*mu_abs,1e-9)) .* max(orf.d_o,1e-9);
            Cd     = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
            dPcalc = 0.5*rho .* ( qmag ./ max(Cd.*Ao,1e-9) ).^2;
            F_lin  = k_sd.*drift + c_lam.*dvel;
            p_up   = orf.p_amb + abs(F_lin)./max(Ap,1e-12);
            dPcav  = max( (p_up - orf.p_cav_eff).*orf.cav_sf, 0 );
            dP_orf = softmin(dPcalc, dPcav, 1e5);          % soft-min clamp
            Q      = Ap * sqrt(dvel.^2 + orf.veps^2);       % piston debisi büyüklüğü
            P_orf  = dP_orf .* Q;
        else
            P_orf = 0*dvel;
        end

        P_sum = sum(P_visc + P_orf, 2);         % toplam güç (t anında)

        % Termal kapasite ve ΔT entegrasyonu (konvektif kayıplı, ileri-Euler)
        nStories = n-1;
        nDtot    = nStories * n_dampers_per_story;
        V_oil_per = resFactor * (Ap * (2*Lgap));
        m_oil_tot   = nDtot * (rho * V_oil_per);
        m_steel_tot = steel_to_oil_mass_ratio * m_oil_tot;
        C_th = max(m_oil_tot*cp_oil + m_steel_tot*cp_steel, eps);

        Nt  = numel(t);
        Tser= zeros(Nt,1); Tser(1)= T0_C;
        dtv = diff(t);
        for k=1:Nt-1
            Pk = 0.5*(P_sum(k) + P_sum(k+1));
            Tser(k+1) = Tser(k) + dtv(k) * ( Pk/C_th - (thermal.hA_W_perK/C_th)*(Tser(k) - thermal.T_env_C) );
            % anlık clamp (süzülmüş taşmaları kes)
            Tser(k+1) = min(max(Tser(k+1), T0_C), T0_C + thermal.dT_max);
        end
        dT_new = Tser(end) - T0_C;
        % son değer clamp
        dT_new = min(max(dT_new, 0), thermal.dT_max);

        % yakınsama ölçütü (stabilite)
        if abs(dT_new - dT_est) <= thermal.tol_K
            dT_est = dT_new;
            break
        end
        % relaxation
        dT_est = thermal.relax*dT_new + (1-thermal.relax)*dT_est;
        T      = T0_C + dT_est;
        % (opsiyonel) bir sonraki iterasyon için başlangıç durumu:
        % z0 = z(end,:).';
    end

    % Son ivme
    a = ( -(M\(C*v.' + K*x.' + dev_force(x,v,k_sd,c_lam,use_orf,orf,rho,Ap,Ao,Qcap,mu_abs))).' ...
          - ag.*r.' );

    % Diag
    diag = struct('dT_est',dT_est,'c_lam',c_lam);

    % ----- iç yardımcılar -----
    function dz = rhs(tt,zz,c_lam_loc,mu_abs_loc)
        x_ = zz(1:n); v_ = zz(n+1:end);
        Fd = dev_force(x_,v_, k_sd,c_lam_loc, use_orf,orf,rho,Ap,Ao,Qcap,mu_abs_loc);
        dv = M \ ( -C*v_ - K*x_ - Fd - M*r*agf(tt) );
        dz = [v_; dv];
    end
end

function F = dev_force(x,v, k_sd,c_lam, use_orf,orf,rho,Ap,Ao,Qcap,mu_abs)
    % x,v: (n x 1) veya (Nt x n) olabilir → çıktı (n x 1) veya (n x Nt)
    singleVec = isvector(x);
    if singleVec, x=x(:).'; v=v(:).'; end
    Nt = size(x,1); n = size(x,2);

    drift = x(:,2:end) - x(:,1:end-1);     % (Nt x n-1)
    dvel  = v(:,2:end) - v(:,1:end-1);

    % Lineer kısım (her iki modda da var)
    F_lin_story = k_sd*drift + c_lam*dvel; % (Nt x n-1)

    if ~use_orf
        F_story = F_lin_story;
    else
        % Akış büyüklüğü (satürasyonlu, yumuşatılmış hız)
        qmag = Qcap * tanh( (Ap/Qcap) * sqrt(dvel.^2 + orf.veps^2) );

        % --- Re düzeltmesi: d_o ile ölçek (μ_abs ile) ---
        Re = (rho .* qmag ./ max(Ao*mu_abs, 1e-9)) .* max(orf.d_o,1e-9);

        % Cd(Re) geçiş eğrisi
        Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);

        % Orifis basınç kaybı
        dP_calc = 0.5*rho .* ( qmag ./ max(Cd.*Ao, 1e-9) ).^2;

        % Kavitasyon üst sınırı (soft-min için tavan)
        p_up   = orf.p_amb + abs(F_lin_story)./max(Ap,1e-12);
        dP_cav = max( (p_up - orf.p_cav_eff) .* orf.cav_sf, 0 );

        % Soft-min klamp (köşeleri yumuşat)
        dP_orf = softmin(dP_calc, dP_cav, 1e5);

        % İşaret
        sgn = dvel ./ sqrt(dvel.^2 + orf.veps^2);

        % Orifis katkısı (hikâye kuvveti)
        F_orf_story = dP_orf .* Ap .* sgn;

        F_story = F_lin_story + F_orf_story;
    end

    % Hikâye → düğümler (eşdeğer yük dağılımı)
    % (alt düğüme +, üst düğüme -; iç düğümler fark)
    F = zeros(Nt,n);
F(:,1)      = -F_story(:,1);
if n>2, F(:,2:n-1) = F_story(:,1:end-1) - F_story(:,2:end); end
F(:,n)      =  F_story(:,end);


    F = F.';                                 % (n x Nt) veya (n x 1)
    if singleVec, F = F(:,1); end
end

function y = softmin(a,b,epsm)
    % C2-sürekli yaklaşık min(a,b)  
    y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
end

function s = tern(c,a,b)
    if c, s=a; else, s=b; end
end

function [t5,t95] = arias_win(t,ag,p1,p2)
    if nargin<3, p1=0.05; p2=0.95; end
    a2 = ag(:).^2; dt = [diff(t); median(diff(t))];
    E  = cumsum(0.5*(a2 + [a2(2:end); a2(end)]).*dt);
    t5  = interp1(E,t,p1*E(end),'linear','extrap');
    t95 = interp1(E,t,p2*E(end),'linear','extrap');
end
