function records = load_ground_motions()
%LOAD_GROUND_MOTIONS Load and preprocess multiple ground-motion records.
%   records = LOAD_GROUND_MOTIONS() loads all ground-motion time histories
%   from <acc_matrix.mat> whose variables follow the pattern acc_matrix1,
%   acc_matrix2, ... The function performs basic checks and preprocessing
%   (demean/detrend and a small high-pass filter) and returns a struct array
%   with time vector, acceleration, and metadata (name, dt, duration,
%   rough PGA/PGV).
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

raw = load('acc_matrix.mat');
fn  = fieldnames(raw);
records = struct('name',{},'t',{},'ag',{},'dt',{},'duration',{},'PGA',{},'PGV',{});

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
                            'PGA', PGA, 'PGV', PGV); %#ok<AGROW>
end

% Print a quick summary table
fprintf('Loaded %d ground-motion records:\n', numel(records));
for k = 1:numel(records)
    r = records(k);
    fprintf('%2d) %-12s dt=%6.4f s dur=%6.2f s PGA=%7.3f PGV=%7.3f\n', ...
        k, r.name, r.dt, r.duration, r.PGA, r.PGV);
end
end
