function [records, scaled] = load_ground_motions(T1, targetIM)
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
    zeta = 0.05;  % default damping ratio
    for k = 1:numel(records)
        records(k).IM = calc_psa(records(k).t, records(k).ag, T1, zeta);
    end

    if nargin < 2 || isempty(targetIM)
        targetIM = mean([records.IM]);
    end

    scaled = records;  % start with a copy
    for k = 1:numel(records)
        s = targetIM / records(k).IM;
        scaled(k).ag   = s * records(k).ag;
        scaled(k).PGA  = s * records(k).PGA;
        scaled(k).PGV  = s * records(k).PGV;
        scaled(k).scale = s;
        scaled(k).IM    = calc_psa(records(k).t, scaled(k).ag, T1, zeta);
    end

    err = abs([scaled.IM] - targetIM) / targetIM * 100;
    fprintf('Target IM = %.3f (PSA @ T1=%.2fs). Max error = %.2f%%\n', ...
        targetIM, T1, max(err));
end
end
