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
use_orifice = true;     % Orifis modeli aç/kapa
use_thermal = true;     % Termal döngü (ΔT ve c_lam(T)) aç/kapa

%% 1–3) Parametrelerin yüklenmesi
% Yapı, damper ve akış/termal parametreleri ayrı bir dosyada tutulur.
parametreler;           % T1 değeri burada hesaplanır

%% 0) Deprem girdisi (ham ivme, m/s^2)
use_scaled = true;               % true → ölçekli kayıt, false → ham kayıt
if use_scaled
    [recs_raw, recs] = load_ground_motions(T1);   % kayıtları yükle ve ölçekle
else
    [recs_raw, ~] = load_ground_motions();        % sadece ham kayıtları yükle
    recs = recs_raw;
end
irec = 1;                        % kullanılacak kayıt indeksi
t    = recs(irec).t;
ag   = recs(irec).ag;

% (Sadece gösterim) Arias %5–%95 penceresi
[t5,t95] = arias_win(t,ag,0.05,0.95);

%% 4) Çözümler
[x0,a0]    = lin_MCK(t,ag,M,C0,K);  % dampersiz

% Lineer (termal kapalı tutarak)
[x_lin,a_lin,diag_lin] = mck_with_damper( ...
    t,ag,M,C0,K, k_sd, c_lam0, false, orf, rho, Ap, A_o, Qcap_big, mu_ref, ...
    false, thermal, T0_C, T_ref_C, b_mu, c_lam_min, c_lam_cap, Lgap, ...
    cp_oil, cp_steel, steel_to_oil_mass_ratio, toggle_gain, story_mask, ...
    n_dampers_per_story, resFactor, cfg);

% Orifisli (+ termal anahtarına göre)
[x_orf,a_orf,diag_orf] = mck_with_damper( ...
    t,ag,M,C0,K, k_sd, c_lam0, use_orifice, orf, rho, Ap, A_o, Qcap_big, mu_ref, ...
    use_thermal, thermal, T0_C, T_ref_C, b_mu, c_lam_min, c_lam_cap, Lgap, ...
    cp_oil, cp_steel, steel_to_oil_mass_ratio, toggle_gain, story_mask, ...
    n_dampers_per_story, resFactor, cfg);

x10_0   = x0(:,10);
x10_lin = x_lin(:,10);
x10_orf = x_orf(:,10);

% 10. kat mutlak ivmeler
a10_0   = a0(:,10)   + ag;
a10_lin = a_lin(:,10)+ ag;
a10_orf = a_orf(:,10)+ ag;

%% Self-check \zeta_{1}
% Toggle ve damper çoğulluğu dikkate alınarak, lineer ve orifis/termal
% senaryoları için birinci mod sönüm oranını hesapla.
nStories = n - 1;
Rvec = toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
mask = story_mask(:); if numel(mask)==1, mask = mask*ones(nStories,1); end
ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
multi = mask .* ndps;                % sütun vektörü

% Damper kaynaklı katkı matrisleri
Kadd = zeros(n);
Cl_add = zeros(n);
Co_add = zeros(n);
for i = 1:nStories
    idx = [i, i+1];
    k_eq  = k_sd * (Rvec(i)^2) * multi(i);
    c_eq_l = diag_lin.c_lam * (Rvec(i)^2) * multi(i);
    c_eq_o = diag_orf.c_lam * (Rvec(i)^2) * multi(i);
    kM = k_eq  * [1 -1; -1 1];
    cM_l = c_eq_l * [1 -1; -1 1];
    cM_o = c_eq_o * [1 -1; -1 1];
    Kadd(idx,idx)  = Kadd(idx,idx)  + kM;
    Cl_add(idx,idx)= Cl_add(idx,idx)+ cM_l;
    Co_add(idx,idx)= Co_add(idx,idx)+ cM_o;
end

K_tot = K + Kadd;
C_lin = C0 + Cl_add;
C_orf = C0 + Co_add;

% Birinci mod için özdeğer/özvektör
[V,D] = eig(K_tot,M);
[w2,ord] = sort(diag(D),'ascend');
phi1 = V(:,ord(1));
w1 = sqrt(w2(1));
normM = phi1.' * M * phi1;
zeta_lin = (phi1.' * C_lin * phi1) / (2*w1*normM);
zeta_orf = (phi1.' * C_orf * phi1) / (2*w1*normM);

%% 5) Grafiklerin çizimi ve kısa özet
grafik;

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
    cp_oil,cp_steel, steel_to_oil_mass_ratio, toggle_gain, story_mask, ...
    n_dampers_per_story, resFactor, cfg)

    n = size(M,1); r = ones(n,1);
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0  = zeros(2*n,1);
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    % --- Toggle ve damper çoğulluğu vektörleri ---
    nStories = n-1;
    Rvec  = toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
    mask  = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps  = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi = (mask .* ndps).';           % satır vektörü
    Rvec  = Rvec.';                     % satır vektörü
    Nvec  = 1:nStories;                 % alt düğümler
    Mvec  = 2:n;                        % üst düğümler

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
        odef = @(tt,z) rhs(tt,z,c_lam,mu_abs,cfg);
        sol  = ode15s(odef,[t(1) t(end)],z0,opts);
        z    = deval(sol,t).';
        x    = z(:,1:n); v = z(:,n+1:end);

        % ---- Termal güncelle (kapalıysa atla) ----
        if ~use_thermal
            break
        end

        % Story farkları
        drift = x(:,Mvec) - x(:,Nvec);      % [Nt x (n-1)]
        dvel  = v(:,Mvec) - v(:,Nvec);

        % Toggle → piston dönüşümü
        drift_p = drift .* Rvec;            % piston yer değiştirmesi
        dvel_p  = dvel  .* Rvec;            % piston hızı

        % Viskoz güç (damper başına)
        P_visc_per = c_lam * (dvel_p.^2);

        % Orifis güç (varsa, damper başına)
        if use_orf
            qmag   = Qcap * tanh( (Ap/Qcap)*sqrt(dvel_p.^2 + orf.veps^2) );
            % Re = (rho*q/Ao/μ)*d_o  → μ_abs ile
            Re     = (rho .* qmag ./ max(Ao*mu_abs,1e-9)) .* max(orf.d_o,1e-9);
            Cd     = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
            dPcalc = 0.5*rho .* ( qmag ./ max(Cd.*Ao,1e-9) ).^2;
            F_lin_p = k_sd.*drift_p + c_lam.*dvel_p;
            p_up   = orf.p_amb + abs(F_lin_p)./max(Ap,1e-12);
            dPcav  = max( (p_up - orf.p_cav_eff).*orf.cav_sf, 0 );
            dP_orf = softmin(dPcalc, dPcav, 1e5);          % soft-min clamp
            Q      = Ap * sqrt(dvel_p.^2 + orf.veps^2);       % piston debisi büyüklüğü
            P_orf_per  = dP_orf .* Q;
        else
            P_orf_per = 0*dvel_p;
        end

        % Toplam güç (t anında, tüm damperler)
        P_sum = sum( (P_visc_per + P_orf_per) .* multi, 2 );

        % Termal kapasite ve ΔT entegrasyonu (konvektif kayıplı, ileri-Euler)
        nStories = n-1;
        nDtot    = sum(multi);
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
    a = ( -(M\(C*v.' + K*x.' + dev_force(t,x,v,k_sd,c_lam,use_orf,orf,rho,Ap,Ao,Qcap,mu_abs,toggle_gain,story_mask,n_dampers_per_story,cfg))).' ...
          - ag.*r.' );
% Diag
    diag = struct('dT_est',dT_est,'c_lam',c_lam);
    
    % ----- iç yardımcılar -----
    function dz = rhs(tt,zz,c_lam_loc,mu_abs_loc,cfg)
        x_ = zz(1:n); v_ = zz(n+1:end);
        Fd = dev_force(tt,x_,v_, k_sd,c_lam_loc, use_orf,orf,rho,Ap,Ao,Qcap,mu_abs_loc,toggle_gain,story_mask,n_dampers_per_story,cfg);
        dv = M \ ( -C*v_ - K*x_ - Fd - M*r*agf(tt) );
        dz = [v_; dv];
    end
end

function F = dev_force(tt,x,v, k_sd,c_lam, use_orf,orf,rho,Ap,Ao,Qcap,mu_abs,toggle_gain,story_mask,n_dampers_per_story,cfg)
    % tt: zaman (skaler veya vektör)
    % x,v: (n x 1) veya (Nt x n) olabilir → çıktı (n x 1) veya (n x Nt)
    singleVec = isvector(x);
    if singleVec, x=x(:).'; v=v(:).'; tt = tt(:); end
    Nt = size(x,1); n = size(x,2);

    nStories = n-1;
    Nvec = 1:nStories; Mvec = 2:n;
    Rvec  = toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
    mask  = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps  = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi = (mask .* ndps).';           % satır vektörü
    Rvec  = Rvec.';                     % satır vektörü

    drift = x(:,Mvec) - x(:,Nvec);     % (Nt x n-1)
    dvel  = v(:,Mvec) - v(:,Nvec);

    % Toggle → piston
    drift_p = drift .* Rvec;
    dvel_p  = dvel  .* Rvec;

    % Lineer kısım piston domaininde (damper başına)
    F_lin_p = k_sd*drift_p + c_lam*dvel_p; % (Nt x n-1)

    if ~use_orf
        F_orf_p = 0*dvel_p;
        F_p = F_lin_p;
    else
        % Akış büyüklüğü (satürasyonlu, yumuşatılmış hız)
        qmag = Qcap * tanh( (Ap/Qcap) * sqrt(dvel_p.^2 + orf.veps^2) );

        % --- Re düzeltmesi: d_o ile ölçek (μ_abs ile) ---
        Re = (rho .* qmag ./ max(Ao*mu_abs, 1e-9)) .* max(orf.d_o,1e-9);

        % Cd(Re) geçiş eğrisi
        Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);

        % Orifis basınç kaybı
        dP_calc = 0.5*rho .* ( qmag ./ max(Cd.*Ao, 1e-9) ).^2;

        % Kavitasyon üst sınırı (soft-min için tavan)
        p_up   = orf.p_amb + abs(F_lin_p)./max(Ap,1e-12);
        dP_cav = max( (p_up - orf.p_cav_eff) .* orf.cav_sf, 0 );

        % Soft-min klamp (köşeleri yumuşat)
        dP_orf = softmin(dP_calc, dP_cav, 1e5);

        % İşaret
        sgn = dvel_p ./ sqrt(dvel_p.^2 + orf.veps^2);

        % Orifis katkısı (piston kuvveti)
        F_orf_p = dP_orf .* Ap .* sgn;

        F_p = F_lin_p + F_orf_p;
    end

    % PF filtresi: sadece rezistif bileşen
    dp_pf = (c_lam*dvel_p + F_orf_p) ./ Ap;
    if isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
        s = tanh(20*dvel_p);
        dp_pf = s .* max(0, s .* dp_pf);
    end
    w_pf = pf_weight(tt, cfg) * cfg.PF.gain;    % Nt x 1
    F_p = k_sd*drift_p + (w_pf .* dp_pf) * Ap;

    % Piston kuvveti → hikâye kuvveti (toggle & çoğulluk)
    F_story = F_p .* (Rvec .* multi);

    % Hikâye → düğümler (eşdeğer yük dağılımı)
    F = zeros(Nt,n);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;

    F = F.';                                 % (n x Nt) veya (n x 1)
    if singleVec, F = F(:,1); end
end

function w = pf_weight(t, cfg)
    w = cfg.on.pressure_force * (1 - exp(-max(t - cfg.PF.t_on, 0) ./ max(cfg.PF.tau, 1e-6)));
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
