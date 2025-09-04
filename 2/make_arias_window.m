function win = make_arias_window(t, ag, varargin)
%MAKE_ARIAS_WINDOW Determine time window capturing main shaking.
%   WIN = MAKE_ARIAS_WINDOW(T, AG) computes the cumulative Arias intensity
%   of the acceleration history AG at times T, identifies the times T5 and
%   T95 when the normalized intensity reaches 5%% and 95%%, and returns a
%   struct describing a window around the significant shaking. Optional
%   parameters can adjust the bounds and padding.
%
%   WIN = MAKE_ARIAS_WINDOW(...,'p1',P1,'p2',P2,'pad',PAD) allows the
%   fractional intensity levels P1/P2 (defaults 0.05/0.95) and the padding
%   duration PAD (default 0.5 s) to be specified. If the shaking duration
%   (T95-T5) is shorter than 5 s the padding is reduced to 0.25 s.
%
%   Output fields of WIN:
%       t5, t95        - times at which cumulative intensity reaches P1,P2
%       pad            - padding actually used [s]
%       t_start,t_end  - start and end times of window [s]
%       idx            - logical index vector for samples within window
%       coverage       - fraction of total Arias intensity in window
%       flag_low_arias - true if coverage < 0.90

p = inputParser;
p.addParameter('p1',0.05,@(x)isscalar(x) && x>=0 && x<=1);
p.addParameter('p2',0.95,@(x)isscalar(x) && x>=0 && x<=1);
p.addParameter('pad',0.5,@(x)isscalar(x) && x>=0);
p.parse(varargin{:});

p1 = p.Results.p1;
p2 = p.Results.p2;
pad = p.Results.pad;

% cumulative Arias intensity (unnormalized)
IA = cumtrapz(t, ag.^2);
IA_tot = IA(end);
IA_norm = IA / IA_tot;

% interpolate times for specified intensity fractions
% ensure monotonic input
[t_unique, iu] = unique(IA_norm);
% corresponding times
interp_t = t(iu);

t5  = interp1(t_unique, interp_t, p1, 'linear');
t95 = interp1(t_unique, interp_t, p2, 'linear');

% adjust padding for short records
dur = t95 - t5;
if dur < 5
    pad = 0.25;
end

t_start = max(t(1),  t5  - pad);
t_end   = min(t(end), t95 + pad);
idx = (t >= t_start) & (t <= t_end);

coverage = trapz(t(idx), ag(idx).^2) / IA_tot;
flag_low_arias = coverage < 0.90;

win = struct('t5',t5,'t95',t95,'pad',pad, ...
             't_start',t_start,'t_end',t_end,'idx',idx, ...
             'coverage',coverage,'flag_low_arias',flag_low_arias);
end
