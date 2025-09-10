%% ================================================================
%  Parametre Tanımları
%  --------------------
%  Bu dosya, ana betik olan <code>damperlinon.m</code> tarafından
%  çağrılır ve analizde kullanılan tüm yapısal, damper ve akış/termal
%  parametrelerini tanımlar. 3 bölümden oluşur:
%    1) 10 katlı çerçevenin kütle, rijitlik ve sönüm özellikleri
%    2) Damper geometrisi ve malzeme özellikleri
%    3) Orifis ve termal modele ilişkin katsayılar
%% ================================================================

%% --- 1) Yapı (10 kat) ---
n  = 10;                               % Kat sayısı
m  = 360e3 * ones(n,1);                % Kat kütleleri [kg]
k  = 6.5e8 * ones(n,1);                % Kat rijitlikleri [N/m]
c  = 6.2e6 * ones(n,1);                % Rayleigh eşdeğeri sönüm [N·s/m]
story_height = 3.0;                    % Kat yüksekliği [m]
% --- Kütle, rijitlik ve sönüm matrisleri ---
M  = diag(m);
K  = zeros(n); C0 = zeros(n);
for i = 1:n
    kL = k(i);    cL = c(i);
    if i < n
        kU = k(i+1); cU = c(i+1);
    else
        kU = 0;    cU = 0;
    end
    K(i,i)   = kL + kU;   C0(i,i)   = cL + cU;
    if i > 1
        K(i,i-1) = -kL;  C0(i,i-1) = -cL;
    end
    if i < n
        K(i,i+1) = -kU;  C0(i,i+1) = -cU;
    end
end

% Birinci doğal periyot (dampersiz referans)
[~,D] = eig(K,M);
w  = sqrt(sort(diag(D),'ascend'));
T1 = 2*pi/w(1);

%% --- 2) Damper geometrisi ve malzeme ---
Dp     = 0.125;     % Piston çapı [m]
Dp_vals   = Dp;     % Tarama aralığı
Lgap   = 0.055;     % Dış gövde/piston aralığı [m]
Lgap_vals = Lgap;   % Tarama aralığı
d_o    = 1.5e-3;    % Orifis çapı [m]
d_o_vals = d_o;     % Tarama aralığı
Lori   = 0.10;      % Orifis uzunluğu [m]
Lori_vals = Lori;   % Tarama aralığı
mu_ref = 0.9;       % Referans viskozite [Pa·s]
mu_ref_vals = mu_ref; % Tarama aralığı

Kd     = 1.6e9;     % Hidrolik sertlik katsayısı [N/m]
Ebody  = 2.1e11;    % Gövde elastisite modülü [Pa]
Gsh    = 79e9;      % Körük kayma modülü [Pa]
d_w    = 12e-3;     % Yay teli çapı [m]
d_w_vals = d_w;     % Tarama aralığı
D_m    = 80e-3;     % Yay orta çapı [m]
D_m_vals = D_m;     % Tarama aralığı
n_turn = 8;         % Yay tur sayısı [-]
n_turn_vals = n_turn; % Tarama aralığı

% Türetilen sabitler (lineer eşdeğer)
Ap    = pi*(Dp.^2)/4;                 % Piston alanı [m^2]
k_h   = Kd*(Ap.^2)./Lgap;             % Hidrolik sertlik [N/m]
k_s   = Ebody*Ap./Lgap;               % Gövde sertliği [N/m]
k_hyd = 1./(1./k_h + 1./k_s);         % Seri bağlanmış eşdeğer
k_p   = Gsh*(d_w.^4)./(8*n_turn.*D_m.^3); % Yay (körük) sertliği [N/m]
k_sd  = k_hyd + k_p;                  % Toplam seri damper sertliği [N/m]
c_lam0 = 12*mu_ref*Lori*(Ap.^2)./(d_o.^4); % Laminer eşdeğer sönüm (T0)

%% --- 3) Orifis ve termal parametreleri ---
rho   = 850;       % Yağ yoğunluğu [kg/m^3]
rho_vals = rho;    % Tarama aralığı
n_orf = 6;         % Kat başına orifis sayısı
n_orf_vals = n_orf; % Tarama aralığı
A_o   = n_orf.*(pi*(d_o.^2)/4);   % Toplam orifis alanı [m^2]

% Orifis katsayıları
orf = struct();
orf.Cd0   = 0.61;                    % Re → 0 limitindeki boşalım katsayısı
orf_Cd0_vals = orf.Cd0;              % Tarama aralığı
orf.CdInf = 0.80;                    % Yüksek Re için boşalım katsayısı
orf_CdInf_vals = orf.CdInf;          % Tarama aralığı
orf.Rec   = 3000;                    % Kritik Reynolds sayısı
orf_Rec_vals = orf.Rec;              % Tarama aralığı
orf.p_exp = 1.1;                     % Geçiş eğrisinin eğimi
orf_p_exp_vals = orf.p_exp;          % Tarama aralığı
orf.p_amb = 1.0e5;                   % Ortam basıncı [Pa]
orf_p_amb_vals = orf.p_amb;          % Tarama aralığı
orf.p_cav_eff = 2.0e3;               % Etkin kavitasyon eşiği [Pa]
orf_p_cav_eff_vals = orf.p_cav_eff;  % Tarama aralığı
orf.cav_sf    = 0.90;                % Kavitasyon emniyet katsayısı
orf_cav_sf_vals = orf.cav_sf;        % Tarama aralığı
orf.d_o   = d_o;                     % Re düzeltmesi için çap [m]
orf.veps  = 0.10;                    % Düşük hız yumuşatma [m/s]

% Akış satürasyonu (sayısal kararlılık için)
Qcap_big = max(orf.CdInf*A_o, 1e-9) * sqrt(2*1.0e9/rho);

% Termal model parametreleri
T0_C      = 25;                      % Başlangıç sıcaklığı [°C]
T_ref_C   = 25;                      % Referans sıcaklık [°C]
b_mu      = -0.013;                  % Viskozite-sıcaklık katsayısı
thermal = struct();
thermal.hA_W_perK = 450;             % Konvektif ısı kaybı katsayısı
thermal_hA_vals = thermal.hA_W_perK; % Tarama aralığı
thermal.T_env_C   = 25;              % Ortam sıcaklığı [°C]
thermal.max_iter  = 3;               % ΔT iterasyon sayısı
thermal.tol_K     = 0.5;             % Yakınsama toleransı [K]
thermal.relax     = 0.5;             % Relaxation katsayısı
thermal.dT_max    = 80;              % Maksimum izin verilen ΔT [K]

% Ek kütle/kapasite verileri
steel_to_oil_mass_ratio = 1.5;       % Çelik/yağ kütle oranı
n_dampers_per_story    = [1];        % Kat başına damper adedi (tarama için dizi)
toggle_gain            = [1.6];      % Toggle kazancı (tarama için dizi)
story_mask             = ones(n-1,1);% Kat maskesi; 1=aktif, 0=damper yok
cp_oil   = 1800;                     % Yağın özgül ısısı [J/(kg·K)]
cp_steel = 500;                      % Çeliğin özgül ısısı [J/(kg·K)]
resFactor = 12;                       % Hacim/kapasite ölçeği

% c_lam(T) sınırları
c_lam_cap      = 2e7;                % Üst sınır (cap)
c_lam_min_frac = 0.05;               % Tabandaki oran
c_lam_min_abs  = 1e5;                % Mutlak taban
c_lam_min      = max(c_lam_min_abs, c_lam_min_frac*c_lam0);

% Basınç-kuvvet filtresi (PF) ayarları
cfg = struct();
cfg.PF.mode      = 'ramp';
cfg.PF.tau       = 0.05;
cfg_PF_tau_vals  = cfg.PF.tau;       % Tarama aralığı
cfg.PF.gain      = 1.0;
cfg_PF_gain_vals = cfg.PF.gain;      % Tarama aralığı
cfg.PF.t_on      = 0;
cfg.PF.auto_t_on = true;
cfg.on.pressure_force     = true;
cfg.on.pf_resistive_only = true;  % sadece rezistif (viskoz+orifis) bileşeni filtrele
