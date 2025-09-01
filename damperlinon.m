%% ================================================================
%  10-Katlı Çerçeve — ODE-only (ham ivme, m/s^2)
%  Free-decay tail (kayıt bitince ag=0) + Orifis + Termo-viskoz geri besleme
%  HIZLI SÜRÜM: tek solver=ode23tb, çift hesap yok, dP istatistik dar
%  + TEŞHİS GRAFİKLERİ: drift/kuvvet, dP_orf(t), T(t), mu(t), kümülatif enerji, kavitasyon payı
% ================================================================

clear; clc; close all;

%% --- Model anahtarları ---
use_orifice = true;     % Orifis modeli aç/kapa
use_thermal = false;     % Termal döngü (ΔT ve c_lam(T)) aç/kapa

%% 0) Deprem girdisi (ham ivme, m/s^2)  → KAYIT BÖLÜMÜ
S       = load('acc_matrix.mat','acc_matrix');   % gerekirse path'i değiştirin
t_raw   = S.acc_matrix(:,1);
ag_raw  = S.acc_matrix(:,2);

% tekilleştir + eş-adımlı yeniden örnekleme (sadece KAYIT bölümü için)
[t_raw,iu] = unique(t_raw,'stable'); 
ag_raw     = ag_raw(iu);
dt         = median(diff(t_raw));
t_rec      = (t_raw(1):dt:t_raw(end)).';              % kayıt zaman vektörü
ag_rec     = interp1(t_raw, ag_raw, t_rec, 'linear'); % kayıt ivmesi

% (Sadece gösterim) Arias %5–%95 penceresi (KAYIT üzerinden)
[t5,t95]   = arias_win(t_rec,ag_rec,0.05,0.95);

%% 1) Yapı (10 kat)
n  = 10;
m  = 2.2e6 * ones(n,1);
k  = 2.95e8 * ones(n,1);
c  = 2.55e6 * ones(n,1);
[M,K,C0] = make_KCM(n,m,k,c);

% T1 (dampersiz referans)
[~,D] = eig(K,M); w = sqrt(sort(diag(D),'ascend')); T1 = 2*pi/w(1);

%% === KAYITTAN SONRA SERBEST SALINIM UZATMASI ===
T_tail_sec = max(20, 12*T1);
t_tail     = (t_rec(end)+dt : dt : t_rec(end)+T_tail_sec).';
t          = [t_rec; t_tail];                     % yeni (uzatılmış) zaman vektörü
ag         = [ag_rec; zeros(size(t_tail))];       % kayıt sonrası ag = 0

%% 2) Damper geometrisi ve malzeme (verdiğin değerler)
Dp=0.125;   Lgap=0.055;  d_o=2.7e-3;  Lori=0.10;  mu_ref=0.9;   % mu_ref [Pa·s]
Kd=1.6e9;   Ebody=2.1e11; Gsh=79e9;   d_w=12e-4;  D_m=80e-3; n_turn=8;

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

orf.Cd0        = 0.61;
orf.CdInf      = 0.80;
orf.Rec        = 3000;
orf.p_exp      = 1.1;
orf.p_amb      = 1.0e5;         % ortam basıncı
orf.p_vap      = 3.0e3;         % yağ buhar basıncı (kavitasyon için)
orf.cav_safety = 0.90;          % emniyet katsayısı
orf.d_o        = d_o;           % Re düzeltmesi için çap
orf.veps       = 0.10;          % [m/s] küçük hız yumuşatma

% Akış satürasyonu (sayısal kararlılık için)
Qcap_big = max(orf.CdInf*A_o, 1e-9) * sqrt(2*1.0e9/rho);

% --- Termal model (ΔT relaxation + clamp + c_lam(T) sınırları) ---
T0_C      = 25;       T_ref_C = 25;   b_mu = -0.013;      % μ(T)=μ_ref*exp(b_mu*(T-Tref))
thermal.hA_W_perK = 150;                                          % konvektif ısı kaybı
thermal.T_env_C   = 25;
thermal.max_iter  = 2;    % hız için 2 iterasyon genelde yeterli
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

%% 4) Çözümler (uzatılmış t, ag üzerinde)
[x0,a0]      = lin_MCK(t,ag,M,C0,K);  % dampersiz

% Lineer (termal kapalı)
[x_lin,a_lin] = mck_with_damper( ...
    t,ag,M,C0,K, k_sd, c_lam0, false, orf, rho, Ap, A_o, Qcap_big, mu_ref, ...
    false, thermal, T0_C, T_ref_C, b_mu, c_lam_min, c_lam_cap, Lgap, ...
    cp_oil, cp_steel, steel_to_oil_mass_ratio, n_dampers_per_story, resFactor);

% Orifisli (+ termal anahtarı)  → 3 çıktı (x, a, diag)
[x_orf,a_orf,diag] = mck_with_damper( ...
    t,ag,M,C0,K, k_sd, c_lam0, use_orifice, orf, rho, Ap, A_o, Qcap_big, mu_ref, ...
    use_thermal, thermal, T0_C, T_ref_C, b_mu, c_lam_min, c_lam_cap, Lgap, ...
    cp_oil, cp_steel, steel_to_oil_mass_ratio, n_dampers_per_story, resFactor);

% 10. kat konum
x10_0   = x0(:,10);
x10_lin = x_lin(:,10);
x10_orf = x_orf(:,10);

% 3. ve 10. kat ivmeler (rölatif). Mutlak için: a_*_abs = a_* + ag;
a3_0   = a0(:,3);    a10_0 = a0(:,10);
a3_lin = a_lin(:,3); a10_lin = a_lin(:,10);
a3_orf = a_orf(:,3); a10_orf = a_orf(:,10);

%% 5) Ana cevap grafikleri
figure('Name','10. Kat yer değiştirme — ODE-only (uzatılmış)','Color','w');
plot(t,x10_0 ,'k','LineWidth',1.4); hold on;
plot(t,x10_lin,'b','LineWidth',1.1);
plot(t,x10_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('x_{10}(t) [m]');
title(sprintf('10-Kat | T1=%.3f s | Arias [%.3f, %.3f] s | Tail=%.1f s',T1,t5,t95,T_tail_sec));
legend('Dampersiz','Lineer','Orifisli','Location','best');

figure('Name','3. Kat ivme — ODE-only (uzatılmış)','Color','w');
plot(t,a3_0 ,'k','LineWidth',1.4); hold on;
plot(t,a3_lin,'b','LineWidth',1.1);
plot(t,a3_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('a_{3}(t) [m/s^2]');
title(sprintf('3. Kat | T1=%.3f s | Arias [%.3f, %.3f] s | Tail=%.1f s',T1,t5,t95,T_tail_sec));
legend('Dampersiz','Lineer','Orifisli','Location','best');

figure('Name','10. Kat ivme — ODE-only (uzatılmış)','Color','w');
plot(t,a10_0 ,'k','LineWidth',1.4); hold on;
plot(t,a10_lin,'b','LineWidth',1.1);
plot(t,a10_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('a_{10}(t) [m/s^2]');
title(sprintf('10. Kat | T1=%.3f s | Arias [%.3f, %.3f] s | Tail=%.1f s',T1,t5,t95,T_tail_sec));
legend('Dampersiz','Lineer','Orifisli','Location','best');

% Kısa özet
fprintf('x10_max  (dampersiz)   = %.4g m\n', max(abs(x10_0)));
fprintf('x10_max  (lineer)      = %.4g m\n', max(abs(x10_lin)));
fprintf('x10_max  (orifisli%s)  = %.4g m\n', tern(use_thermal,'+termal',''), max(abs(x10_orf)));
if use_orifice && isfield(diag,'dT_est')
    fprintf('Termal: ΔT_est=%.2f K | c_lam(final)=%.3e N·s/m | dP_orf q50=%.3e q95=%.3e q99=%.3e Pa\n', ...
        diag.dT_est, diag.c_lam, diag.dP_q50, diag.dP_q95, diag.dP_q99);
    fprintf('Kavitasyon payı (zaman*kat): %.1f %%\n', 100*diag.cav_frac_all);
end
fprintf('Not: ag uzatmasıyla serbest salınım %.1f s sürdürülmüştür.\n', T_tail_sec);

%% 6) TEŞHİS GRAFİKLERİ (orifis+termal)
if use_orifice
    % 6a) Hikâye driftleri (tüm katlar)
    figure('Name','Hikâye Driftleri (tüm katlar)','Color','w');
    plot(t, diag.drift_t, 'LineWidth',0.9); grid on;
    xlabel('t [s]'); ylabel('\Delta x_i(t) [m]');
    title('Hikâye Driftleri (\Delta x_i = x_{i+1}-x_i)');
    legend(arrayfun(@(i)sprintf('Kat %d-%d',i,i+1),1:n-1,'uni',0),'Location','eastoutside');

    % 6b) Hikâye damper kuvvetleri (lineer+orifis toplamı)
    figure('Name','Hikâye Damper Kuvvetleri (tüm katlar)','Color','w');
    plot(t, diag.F_story_t, 'LineWidth',0.9); grid on;
    xlabel('t [s]'); ylabel('F_{story,i}(t) [N]');
    title('Hikâye Damper Kuvvetleri (k_{sd}\Delta x + c_{lam}\Delta v + F_{orf})');
    legend(arrayfun(@(i)sprintf('Kat %d-%d',i,i+1),1:n-1,'uni',0),'Location','eastoutside');

    % 6c) dP_orf(t) — zarf (max) ve medyan + q-quantile çizgileri
    figure('Name','Orifis Basınç Düşümü — zaman zarfı','Color','w');
    plot(t, diag.dP_orf_max_t, 'LineWidth',1.2); hold on;
    plot(t, diag.dP_orf_med_t, '--', 'LineWidth',1.0);
    yline(diag.dP_q50,'k--','q50','HandleVisibility','off');
    yline(diag.dP_q95,'r--','q95','HandleVisibility','off');
    yline(diag.dP_q99,'m--','q99','HandleVisibility','off');
    grid on; xlabel('t [s]'); ylabel('\Delta P_{\it orf}(t) [Pa]');
    title('\Delta P_{orf}(t): maks-zarf ve medyan (tüm katlar)');
    legend('max(zarf)','medyan','Location','best');

    % 6d) Kavitasyon payı (zaman içinde, katlar ort.)
    figure('Name','Kavitasyon payı (zaman serisi)','Color','w');
    plot(t, diag.cav_frac_t, 'LineWidth',1.2); grid on;
    ylim([0 1]); xlabel('t [s]'); ylabel('kavitasyon payı (0–1)');
    title('Kavitasyon Klamplama Payı (zaman boyunca, katlar ort.)');
end

if use_thermal
    % 6e) T(t) ve μ(t)
    figure('Name','Termal: T(t) ve \mu(t)','Color','w');
    yyaxis left;
    plot(t, diag.T_series, 'LineWidth',1.2); ylabel('T(t) [^\circC]');
    yyaxis right;
    plot(t, diag.mu_series, '--', 'LineWidth',1.2); ylabel('\mu(t) [Pa·s]');
    grid on; xlabel('t [s]');
    title('Termo-viskoz durum: Sıcaklık ve Viskozite');

    % 6f) Kümülatif enerji
    figure('Name','Kümülatif Enerji (viskoz + orifis + toplam)','Color','w');
    plot(t, diag.E_visc_cum, 'LineWidth',1.2); hold on;
    plot(t, diag.E_orf_cum , 'LineWidth',1.2);
    plot(t, diag.E_tot_cum , 'k', 'LineWidth',1.4);
    grid on; xlabel('t [s]'); ylabel('E(t) [J] (ölçekli)');
    title('Damper Enerjisi — Kümülatif');
    legend('E_{visc}','E_{orf}','E_{toplam}','Location','best');
end

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
    Tspan = t(end)-t(1);
    odef = @(tau,z) Tspan * [ z(n+1:end);
                              -(M\(C*z(n+1:end) + K*z(1:n))) - r*agf(t(1)+tau*Tspan) ];
    z0 = zeros(2*n,1);
    opts = odeset('RelTol',2e-2,'AbsTol',[1e-6*ones(n,1); 1e-5*ones(n,1)], ...
                  'MaxStep',5*median(diff(t))/Tspan,'InitialStep',median(diff(t))/Tspan);
    sol = ode23tb(odef,[0 1],z0,opts);
    z = deval(sol,(t-t(1))/Tspan).';
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

    % --- Termal arayıncı başlangıç ---
    T      = T0_C;     % anlık sıcaklık tahmini
    dT_est = 0;        % toplam artış tahmini

    % Döngü sayısı
    nIter = tern(use_thermal, thermal.max_iter, 1);

    % Tek solver: ode23tb
    Tspan = t(end)-t(1);
    dtm   = median(diff(t));
    opts  = odeset('RelTol',2e-2, ...
                   'AbsTol',[1e-6*ones(n,1); 1e-5*ones(n,1)], ...
                   'InitialStep', max(0.5*dtm/Tspan,1e-6), ...
                   'MaxStep',     5*dtm/Tspan);

    % --- Isı kapasitesi (sabit) ---
    nStories   = n-1;
    nDtot      = nStories * n_dampers_per_story;
    V_oil_per  = resFactor * (Ap * (2*Lgap));
    m_oil_tot  = nDtot * (rho * V_oil_per);
    m_steel_tot= steel_to_oil_mass_ratio * m_oil_tot;
    C_th       = max(m_oil_tot*cp_oil + m_steel_tot*cp_steel, eps);

    % === DIAG EMNİYET ÖN-TAHSİS (boş bile olsa alanlar var) ===
    Nt = numel(t); Zt = zeros(Nt, max(n-1,1)); z1 = zeros(Nt,1);
    diag = struct( ...
        'drift_t',      Zt, ...
        'dvel_t',       Zt, ...
        'F_story_t',    Zt, ...
        'F_lin_story',  Zt, ...
        'dP_orf_time',  Zt, ...
        'dP_orf_max_t', z1, ...
        'dP_orf_med_t', z1, ...
        'P_visc_t',     Zt, ...
        'P_orf_t',      Zt, ...
        'E_visc_cum',   z1, ...
        'E_orf_cum',    z1, ...
        'E_tot_cum',    z1, ...
        'cav_frac_t',   z1, ...
        'cav_frac_all', 0, ...
        'T_series',     T0_C*ones(Nt,1), ...
        'mu_series',    mu_ref*exp(b_mu*(T0_C-T_ref_C))*ones(Nt,1), ...
        'dT_est',       0, ...
        'c_lam',        c_lam0, ...
        'dP_q50',       NaN, ...
        'dP_q95',       NaN, ...
        'dP_q99',       NaN );

    % ---- iteratif μ–T döngüsü ----
    for it = 1:nIter
        % viskozite/sönüm güncelle
        mu_scale = max(1e-6, exp(b_mu*(T - T_ref_C)));
        mu_abs   = mu_ref * mu_scale;             % [Pa·s]
        c_raw    = c_lam0 * mu_scale;
        c_lam    = min(max(c_raw, c_lam_min), c_lam_cap);

        % ODE çözümü
        odef = @(tau,z) rhs_tau(tau,z,c_lam,mu_abs);
        sol  = ode23tb(odef,[0 1],z0,opts);
        z    = deval(sol,(t-t(1))/Tspan).';
        x    = z(:,1:n); v = z(:,n+1:end);

        % Story farkları
        drift = x(:,2:end) - x(:,1:end-1);      % [Nt x (n-1)]
        dvel  = v(:,2:end) - v(:,1:end-1);
        F_lin_story = k_sd.*drift + c_lam.*dvel;

        % Orifis (isteğe bağlı)
        if use_orf
            qmag   = Qcap * tanh( (Ap/Qcap)*sqrt(dvel.^2 + orf.veps^2) );
            Re     = (rho .* qmag ./ max(Ao*mu_abs,1e-9)) .* max(orf.d_o,1e-9);
            Cd     = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
            dPcalc = 0.5*rho .* ( qmag ./ max(Cd.*Ao,1e-9) ).^2;
            p_up   = orf.p_amb + abs(F_lin_story)./max(Ap,1e-12);
            dPcav  = max( (p_up - orf.p_vap).*orf.cav_safety, 0 );
            dP_orf = softmin(dPcalc, dPcav, 1e5);
            sgn    = dvel ./ sqrt(dvel.^2 + orf.veps^2);
            F_orf_story = dP_orf .* Ap .* sgn;
            Qs     = Ap * sqrt(dvel.^2 + orf.veps^2);
            P_orf  = dP_orf .* Qs;
            % kavitasyon payları
            cav_mask     = dPcav < dPcalc;
            cav_frac_t   = mean(cav_mask, 2);
            cav_frac_all = mean(cav_mask(:));
            % dP yüzdelleri (tüm zaman*kat)
            dq           = prctile(dP_orf(:), [50 95 99]);
            dP_q50=dq(1); dP_q95=dq(2); dP_q99=dq(3);
            dP_orf_max_t = max(dP_orf, [], 2);
            dP_orf_med_t = median(dP_orf, 2);
        else
            F_orf_story  = 0*dvel;
            P_orf        = 0*dvel;
            cav_frac_t   = zeros(Nt,1);
            cav_frac_all = 0;
            dP_q50=0; dP_q95=0; dP_q99=0;
            dP_orf       = 0*dvel;
            dP_orf_max_t = zeros(Nt,1);
            dP_orf_med_t = zeros(Nt,1);
        end

        % Toplam hikâye kuvveti ve güçler
        F_story = F_lin_story + F_orf_story;
        P_visc  = c_lam * (dvel.^2);

        % Termal entegrasyon (enerji → sıcaklık)
        P_sum = sum(P_visc + P_orf, 2);
        Tser  = zeros(Nt,1);  Tser(1)= T0_C;
        dtv   = diff(t);
        for k=1:Nt-1
            Pk = 0.5*(P_sum(k) + P_sum(k+1));
            Tser(k+1) = Tser(k) + dtv(k) * ( Pk/C_th - (thermal.hA_W_perK/C_th)*(Tser(k) - thermal.T_env_C) );
            Tser(k+1) = min(max(Tser(k+1), T0_C), T0_C + thermal.dT_max);
        end
        dT_new = min(max(Tser(end) - T0_C, 0), thermal.dT_max);

        % Kümülatif enerjiler
        E_visc_cum = cumtrapz(t, sum(P_visc,2));
        E_orf_cum  = cumtrapz(t, sum(P_orf ,2));
        E_tot_cum  = E_visc_cum + E_orf_cum;

        % ==== HER İTERASYONDA DIAG GÜNCELLE (garanti) ====
        diag.drift_t      = drift;
        diag.dvel_t       = dvel;
        diag.F_story_t    = F_story;
        diag.F_lin_story  = F_lin_story;
        diag.dP_orf_time  = dP_orf;
        diag.dP_orf_max_t = dP_orf_max_t;
        diag.dP_orf_med_t = dP_orf_med_t;
        diag.P_visc_t     = P_visc;
        diag.P_orf_t      = P_orf;
        diag.E_visc_cum   = E_visc_cum;
        diag.E_orf_cum    = E_orf_cum;
        diag.E_tot_cum    = E_tot_cum;
        diag.cav_frac_t   = cav_frac_t;
        diag.cav_frac_all = cav_frac_all;
        diag.T_series     = Tser;
        diag.mu_series    = mu_ref * exp(b_mu*(Tser - T_ref_C));
        diag.dT_est       = dT_new;
        diag.c_lam        = c_lam;
        diag.dP_q50       = dP_q50;
        diag.dP_q95       = dP_q95;
        diag.dP_q99       = dP_q99;

        % Yakınsama testi
        if ~use_thermal || abs(dT_new - dT_est) <= thermal.tol_K
            dT_est = dT_new;
            break
        end
        dT_est = thermal.relax*dT_new + (1-thermal.relax)*dT_est;
        T      = T0_C + dT_est;
        % z0 = z(end,:).';  % istersen sıcak başlat
    end

    % Son ivme (tek kuvvet hesabı)
    a = ( -(M\(C*v.' + K*x.' + dev_force(x,v,k_sd,c_lam,use_orf,orf,rho,Ap,Ao,Qcap,mu_abs))).' ...
          - ag.*r.' );

    % ----- iç yardımcılar -----
    function dz = rhs_tau(tau,zz,c_lam_loc,mu_abs_loc)
        x_ = zz(1:n); v_ = zz(n+1:end);
        Fd = dev_force(x_,v_, k_sd,c_lam_loc, use_orf,orf,rho,Ap,Ao,Qcap,mu_abs_loc);
        dv = M \ ( -C*v_ - K*x_ - Fd - M*r*agf(t(1)+tau*Tspan) );
        dz = Tspan * [v_; dv];
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

        % --- Re düzeltmesi: d_o ile ölçekle ---
        Re = (rho .* qmag ./ max(Ao*mu_abs, 1e-9)) .* max(orf.d_o,1e-9);

        % Cd(Re) geçiş eğrisi
        Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);

        % Orifis basınç kaybı
        dP_calc = 0.5*rho .* ( qmag ./ max(Cd.*Ao, 1e-9) ).^2;

        % Kavitasyon üst sınırı (soft-min için tavan)
        p_up   = orf.p_amb + abs(F_lin_story)./max(Ap,1e-12);
        dP_cav = max( (p_up - orf.p_vap) .* orf.cav_safety, 0 );

        % Soft-min klamp (köşeleri yumuşat)
        dP_orf = softmin(dP_calc, dP_cav, 1e5);

        % İşaret
        sgn = dvel ./ sqrt(dvel.^2 + orf.veps^2);

        % Orifis katkısı (hikâye kuvveti)
        F_orf_story = dP_orf .* Ap .* sgn;

        F_story = F_lin_story + F_orf_story;
    end

    % Hikâye → düğümler (G^T çarpımı eşdeğeri)
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
function s = tern(c,a,b), if c, s=a; else, s=b; end, end

function [t5,t95] = arias_win(t,ag,p1,p2)
    if nargin<3, p1=0.05; p2=0.95; end
    a2 = ag(:).^2; dt = [diff(t); median(diff(t))];
    E  = cumsum(0.5*(a2 + [a2(2:end); a2(end)]).*dt);
    t5  = interp1(E,t,p1*E(end),'linear','extrap');
    t95 = interp1(E,t,p2*E(end),'linear','extrap');
end
