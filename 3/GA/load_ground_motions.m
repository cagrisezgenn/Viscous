function [records, scaled, meta] = load_ground_motions(T1, opts)
%LOAD_GROUND_MOTIONS Zemin hareketi kayitlarini yukler, on isler ve olcekler.
%   [RAW, SCALED, META] = LOAD_GROUND_MOTIONS(T1, OPTS) fonksiyonu
%   acc_matrix.mat dosyasindaki acc_matrix1, acc_matrix2, ... degiskenlerini
%   okuyarak her kayit icin zaman vektoru, ivme ve temel buyuklukleri
%   iceren bir yapi dizisi dondurur.
%
%   T1 girilirse kayitlarin 5%% sonumlu T1 periyodundaki veya tanimli
%   bir periyot bandindaki yapay spektral ivmeleri (PSA) hesaplanir ve
%   tum kayitlar ortak bir hedef IM degerine olceklendirilir.
%
%   OPTS alanlari:
%       hp_cut   - yuksek gecis filtresi frekansi [Hz]
%       IM_mode  - 'T1' veya 'band'
%       band_fac - T1 cevresindeki band carpani [alt ust]
%       band_N   - bant icindeki periyot sayisi
%       s_bounds - olcek katsayisi sinirlari [min max]
%
%   Ciktilar:
%       records - ham kayitlar
%       scaled  - olceklenmis kayitlar
%       meta    - islemlere dair ozet bilgiler

%% Girdi Parametreleri
if nargin < 2, opts = struct(); end
hp_cut   = Utils.getfield_default(opts,'hp_cut',0.05);   % yuksek gecis [Hz]
IM_mode  = Utils.getfield_default(opts,'IM_mode','band');
band_fac = Utils.getfield_default(opts,'band_fac',[0.8 1.2]);
band_N   = Utils.getfield_default(opts,'band_N',15);
s_bounds = Utils.getfield_default(opts,'s_bounds',[0.3 1.2]);

%% MAT dosyasını Yükle
raw = load('acc_matrix.mat');
fn  = fieldnames(raw);
records = struct('name',{},'t',{},'ag',{},'dt',{},'duration',{},'PGA',{},'PGV',{},'IM',{},'scale',{});

%% Ön İşleme
for k = 1:numel(fn)
    A = raw.(fn{k});
    t  = A(:,1);   ag = A(:,2);

    % 3A) Zaman vektorunu duzenle
    [t,iu] = unique(t,'stable'); ag = ag(iu);
    dt = median(diff(t));
    assert(max(abs(diff(t) - dt)) < 1e-6, 'Zaman örnekleme aralığı düzensiz');
    t  = (t(1):dt:t(end)).';
    ag = interp1(A(:,1), A(:,2), t, 'linear');

    % 3B) Temel birim ve filtreleme kontrolleri
    assert(max(abs(ag)) < 100, 'İvme büyüklüğü birim hatasına işaret ediyor');
    ag = detrend(ag,0);             % ortalama
    ag = detrend(ag,1);             % dogrusal trend
    try
        Fs = 1/dt;
        Wn = hp_cut/(Fs/2);
        [b,a] = butter(2,Wn,'high');
        ag = filtfilt(b,a,ag);
    catch
        % Butter/Filtfilt yoksa yalnizca detrend uygulanir
    end

    % 3C) Metaveri
    duration = t(end) - t(1);
    v  = cumtrapz(t,ag);
    PGA = max(abs(ag));
    PGV = max(abs(v));

    records(end+1) = struct('name',fn{k},'t',t,'ag',ag,'dt',dt, ...
                             'duration',duration,'PGA',PGA,'PGV',PGV, ...
                             'IM',[],'scale',1); %#ok<AGROW>
end

%% Yükleme Özeti
fprintf('Toplam %d zemin hareketi kaydı yüklendi:\n', numel(records));
for k = 1:numel(records)
    r = records(k);
    fprintf('%2d) %-12s dt=%6.4f s dur=%6.2f s PGA=%7.3f PGV=%7.3f\n', ...
        k, r.name, r.dt, r.duration, r.PGA, r.PGV);
end

scaled = [];
meta = struct();

%% IM Hesabı ve Ölçekleme (T1 sağlandığında)
if nargin >= 1 && ~isempty(T1)
    %% IM Hesabı
    for k = 1:numel(records)
        records(k).IM = compute_IM(records(k).t, records(k).ag, IM_mode, T1, band_fac, band_N);
    end

    %% Hedef IM Seçimi ve TRIM
    IM = [records.IM];
    IM_low  = max(s_bounds(1)*IM);
    IM_high = min(s_bounds(2)*IM);
    targetIM0 = median(IM);
    targetIM  = min(max(targetIM0, IM_low), IM_high);
    doClip = (IM_low > IM_high);

    max_iter = 3; dropped = {};
    for it = 1:max_iter
        IM_low  = max(s_bounds(1)*IM);
        IM_high = min(s_bounds(2)*IM);
        if IM_low <= IM_high, break; end
        [~,idx] = max(abs(log(IM) - median(log(IM))));
        dropped{end+1} = records(idx).name; %#ok<AGROW>
        records(idx) = [];
        IM(idx) = [];
    end
    if ~isempty(dropped)
        fprintf('TRIM: ayıklanan uç değerler = %s\n', strjoin(dropped,', '));
    end
    IM_low  = max(s_bounds(1)*IM);
    IM_high = min(s_bounds(2)*IM);
    targetIM0 = median(IM);
    targetIM  = min(max(targetIM0, IM_low), IM_high);
    doClip = (IM_low > IM_high);

    %% Ölçekleme
    scaled = records; n_clipped = 0; s_all = zeros(1,numel(records));
    for k = 1:numel(records)
        s_raw = targetIM / records(k).IM;
        if doClip
            s = min(max(s_raw, s_bounds(1)), s_bounds(2));
            if abs(s - s_raw) > 1e-12, n_clipped = n_clipped + 1; end
        else
            s = s_raw;
        end
        s_all(k) = s;

        scaled(k).ag    = s * records(k).ag;
        scaled(k).PGA   = s * records(k).PGA;
        scaled(k).PGV   = s * records(k).PGV;
        scaled(k).scale = s;
        scaled(k).s_clipped = doClip && (abs(s - s_raw) > 1e-12);
        scaled(k).trimmed   = false;

        scaled(k).IM = compute_IM(scaled(k).t, scaled(k).ag, IM_mode, T1, band_fac, band_N);
    end

    %% Hata ve Log
    err = abs([scaled.IM] - targetIM) / max(targetIM, eps) * 100;
    if strcmpi(IM_mode,'band')
        modeStr = 'band';
    else
        modeStr = 'PSA@T1';
    end
    fprintf('Hedef IM = %.3f (%s). Maks hata = %.2f%% | uygun aralık=[%.3f, %.3f] | s_min=%.2f s_max=%.2f', ...
        targetIM, modeStr, max(err), IM_low, IM_high, min(s_all), max(s_all));
    if doClip
        fprintf(' | KIRPILAN=%d', n_clipped);
    else
        fprintf(' | KIRPILAN=0');
    end
    fprintf('\n');

    meta.IM_mode  = IM_mode;
    meta.band_fac = band_fac;
    meta.s_bounds = s_bounds;
    if exist('dropped','var'), meta.TRIM_names = dropped; else, meta.TRIM_names = {}; end
end
end

%% ==== Yerel Fonksiyonlar ====
function IM = compute_IM(t, ag, mode, T1, band_fac, band_N)
%COMPUTE_IM Kayit icin hedef IM degerini hesaplar.
%   IM = COMPUTE_IM(T, AG, MODE, T1, BAND_FAC, BAND_N) fonksiyonu,
%   verilen zaman vektoru t ve ivme ag icin, secilen modda
%   yapay spektral ivme (IM) dondurur.

zeta = 0.05;
if strcmpi(mode,'band')
    Tgrid = linspace(band_fac(1)*T1, band_fac(2)*T1, band_N);
    Sa = zeros(size(Tgrid));
    for i = 1:numel(Tgrid)
        Sa(i) = calc_psa(t, ag, Tgrid(i), zeta);
    end
    IM = exp(mean(log(Sa + eps)));
else
    IM = calc_psa(t, ag, T1, zeta);
end
end

function Sa = calc_psa(t, ag, T, zeta)
%CALC_PSA Tek bir kayit icin yapay spektral ivme hesabi.
%   Sa = CALC_PSA(t, ag, T, zeta) fonksiyonu, T periyotlu ve zeta
%   sonumlu tek serbestlik dereceli osilatorun mutlak ivme tepkisinin
%   maksimum degerini dondurur.

w = 2*pi / T;
agf = griddedInterpolant(t, ag, 'linear', 'nearest');

odef = @(tt, y)[ y(2);
                 -2*zeta*w*y(2) - w*w*y(1) - agf(tt) ];

y0 = [0;0];
[~, y] = ode45(odef, t, y0);

x  = y(:,1);
xd = y(:,2);
xdd = -2*zeta*w*xd - w*w*x - ag;
abs_acc = xdd + ag;
Sa = max(abs(abs_acc));
end
