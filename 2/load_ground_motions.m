function [records, scaled] = load_ground_motions(T1, ~)
%LOAD_GROUND_MOTIONS Load and preprocess multiple ground-motion records.
%   R = LOAD_GROUND_MOTIONS() loads all ground-motion time histories from
%   <acc_matrix.mat> whose variables follow the pattern acc_matrix1,
%   acc_matrix2, ... The function performs basic checks and preprocessing
%   (demean/detrend and a small high-pass filter) and returns a struct
%   array with time vector, acceleration, and metadata (name, dt, duration,
%   rough PGA/PGV).
%
%   [RAW, SCALED] = LOAD_GROUND_MOTIONS(T1, TARGETIM) additionally computes
%   an intensity measure (PSA at period T1 with 5%% damping) for each record
%   and scales the accelerations so that the spectral acceleration equals
%   TARGETIM. The unscaled set is returned as RAW and the scaled set as
%   SCALED. The scaling factor is stored in field <scale>.
%   If TARGETIM is empty, the mean IM of the raw set is used.
%
%   Requirements: Signal Processing Toolbox for BUTTER/FILTFILT. If not
%   available, preprocessing will fall back to mean removal only.
%
%   Output fields for each record:
%       name     - variable name inside the MAT-file
%       t, ag    - time vector and acceleration (m/s^2)
%       dt       - sampling interval (s)
%       duration - total duration (s)
%       PGA      - peak ground acceleration (m/s^2)
%       PGV      - peak ground velocity (m/s)
%       IM       - (optional) PSA at T1, 5%% damping [m/s^2]
%       scale    - (scaled set only) scale factor s such that IM*s = TARGETIM

hp_cut = 0.05;      % high-pass [Hz] (var)
IM_mode = 'band';     % 'T1'  ya da 'band'  (band=0.8–1.2*T1 ortalaması)
band_fac = [0.8 1.2];
band_N   = 15;      % band içi T sayısı (log-uzay veya lin fark etmez)
s_bounds = [0.5 2];  % ölçek sınırı: s ∈ [0.5, 2.0]

raw = load('acc_matrix.mat');
fn  = fieldnames(raw);
records = struct('name',{},'t',{},'ag',{},'dt',{},'duration',{},'PGA',{},'PGV',{},'IM',{},'scale',{});

hp_cut = 0.05;   % small high-pass corner [Hz]

for k = 1:numel(fn)
    A = raw.(fn{k});
    t  = A(:,1);   ag = A(:,2);

    % Ensure monotonic time and uniqueness
    [t,iu] = unique(t,'stable'); ag = ag(iu);
    dt = median(diff(t));
    % Check for uniform sampling
    assert(max(abs(diff(t) - dt)) < 1e-6, 'Non-uniform sampling interval');
    t  = (t(1):dt:t(end)).';
    ag = interp1(A(:,1), A(:,2), t, 'linear');

    % Basic unit check (expecting m/s^2)
    assert(max(abs(ag)) < 100, 'Acceleration magnitude suggests wrong units');

    % Preprocess: demean/detrend
    ag = detrend(ag, 0);          % remove mean
    ag = detrend(ag, 1);          % remove linear trend

    % High-pass filter (~0.05 Hz). Fall back gracefully if toolbox missing.
    try
        Fs = 1/dt;
        Wn = hp_cut / (Fs/2);
        [b,a] = butter(2, Wn, 'high');
        ag = filtfilt(b,a,ag);
    catch
        % Butter/filtfilt not available; already demeaned/detrended
    end

    % Metadata
    duration = t(end) - t(1);
    v  = cumtrapz(t, ag);  % integrate acceleration
    PGA = max(abs(ag));
    PGV = max(abs(v));

    records(end+1) = struct('name', fn{k}, 't', t, 'ag', ag, ...
                            'dt', dt, 'duration', duration, ...
                            'PGA', PGA, 'PGV', PGV, ...
                            'IM', [], 'scale', 1); %#ok<AGROW>
end

% Print a quick summary table
fprintf('Loaded %d ground-motion records:\n', numel(records));
for k = 1:numel(records)
    r = records(k);
    fprintf('%2d) %-12s dt=%6.4f s dur=%6.2f s PGA=%7.3f PGV=%7.3f\n', ...
        k, r.name, r.dt, r.duration, r.PGA, r.PGV);
end

scaled = [];
if nargin >= 1 && ~isempty(T1)
    zeta = 0.05;

    % --- 2A) Her kayıt için IM hesapla (T1 veya band) ---
    for k = 1:numel(records)
        t  = records(k).t;  ag = records(k).ag;
        if strcmpi(IM_mode,'band')
            % 0.8–1.2*T1 band-ortalama PSA (geometrik ortalama)
            Tgrid = linspace(band_fac(1)*T1, band_fac(2)*T1, band_N);
            Sa = zeros(size(Tgrid));
            for i = 1:numel(Tgrid)
                Sa(i) = calc_psa(t, ag, Tgrid(i), zeta);
            end
            records(k).IM = exp(mean(log(Sa+eps)));  % geo-mean
        else
            % Klasik: PSA @ T1
            records(k).IM = calc_psa(t, ag, T1, zeta);
        end
    end

    % --- 2B) Scale-cap uyumlu hedef IM seçimi ---
    IM = [records.IM];
    % her kayıt için s ∈ [0.5, 2.0] ⇒ targetIM ∈ [0.5*IM_k, 2.0*IM_k]
    IM_low  = max(s_bounds(1)*IM);   % yapılabilir alt sınır
    IM_high = min(s_bounds(2)*IM);   % yapılabilir üst sınır

    % başlangıç hedef: median(IM); sonra yapılabilir aralığa sıkıştır
    targetIM0 = median(IM);
    targetIM  = min(max(targetIM0, IM_low), IM_high);

    % Eğer aralık boş ise (IM_low > IM_high) → son çare: clipping
    doClip = (IM_low > IM_high);
    IM = [records.IM];
max_iter = 3; dropped = {};
for it = 1:max_iter
    IM_low  = max(s_bounds(1)*IM);
    IM_high = min(s_bounds(2)*IM);
    if IM_low <= IM_high
        break
    end
    % Uç kaydı tespit: median(IM)'den en uzak olan
    [~,idx] = max(abs(log(IM) - median(log(IM))));
    dropped{end+1} = records(idx).name; %#ok<AGROW>
    records(idx) = [];   % kaydı at
    IM(idx) = [];
end
if ~isempty(dropped)
    fprintf('TRIM: dropped outliers = %s\n', strjoin(dropped,', '));
end
% Sonra targetIM = min(max(median(IM), IM_low), IM_high) ile devam...


    % --- 2C) Ölçekle ve (opsiyonel) clip’i raporla ---
    scaled = records;  n_clipped = 0;  s_all = zeros(1,numel(records));
    for k = 1:numel(records)
        s_raw = targetIM / records(k).IM;
        if doClip
            s = min(max(s_raw, s_bounds(1)), s_bounds(2)); % clip
            if abs(s - s_raw) > 1e-12, n_clipped = n_clipped + 1; end
        else
            s = s_raw;  % zaten feasible aralıkta
        end
        s_all(k) = s;

        % ölçeklenmiş kayıt
        scaled(k).ag    = s * records(k).ag;
        scaled(k).PGA   = s * records(k).PGA;
        scaled(k).PGV   = s * records(k).PGV;
        scaled(k).scale = s;

        % yeniden IM (T1 veya band ile tutarlılık)
        if strcmpi(IM_mode,'band')
            Tgrid = linspace(band_fac(1)*T1, band_fac(2)*T1, band_N);
            Sa = zeros(size(Tgrid));
            for i = 1:numel(Tgrid), Sa(i) = calc_psa(scaled(k).t, scaled(k).ag, Tgrid(i), zeta); end
            scaled(k).IM = exp(mean(log(Sa+eps)));
        else
            scaled(k).IM = calc_psa(scaled(k).t, scaled(k).ag, T1, zeta);
        end
    end

    % Hata ve log
   % Hata ve log
err = abs([scaled.IM] - targetIM) / max(targetIM, eps) * 100;

% IM modu yazısı (tern yok; if/else kullan)
if strcmpi(IM_mode,'band')
    modeStr = 'band';
else
    modeStr = 'PSA@T1';
end

fprintf('Target IM = %.3f (%s). Max error = %.2f%% | feasible=[%.3f, %.3f] | s_min=%.2f s_max=%.2f', ...
    targetIM, modeStr, max(err), IM_low, IM_high, min(s_all), max(s_all));

if doClip
    fprintf(' | CLIPPED=%d', n_clipped);
end
fprintf('\n');


end
