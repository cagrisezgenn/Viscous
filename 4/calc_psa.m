function Sa = calc_psa(t, ag, T, zeta)
%CALC_PSA Compute pseudo spectral acceleration for a single record.
%   Sa = CALC_PSA(t, ag, T, zeta) returns the peak absolute acceleration
%   response (pseudo spectral acceleration) of a linear single-degree-of-
%   freedom oscillator with period T and damping ratio zeta when subjected
%   to ground acceleration ag defined at times t.
%
%   t    - time vector [s]
%   ag   - ground acceleration vector [m/s^2]
%   T    - oscillator natural period [s]
%   zeta - damping ratio (e.g., 0.05 for 5%%)
%
%   The equation of motion solved is:
%       x'' + 2*zeta*w*x' + w^2*x = -ag(t)
%   where w = 2*pi/T. The absolute acceleration response is x'' + ag.
%
%   Note: Requires base functions (no special toolboxes).

w = 2*pi / T;
% Interpolant for ground acceleration
agf = griddedInterpolant(t, ag, 'linear', 'nearest');

% ODE state: [x; xdot]
odef = @(tt, y)[ y(2);
                 -2*zeta*w*y(2) - w*w*y(1) - agf(tt) ];

% Integrate using ODE45 at specified time points
y0 = [0; 0];
[~, y] = ode45(odef, t, y0);

x = y(:,1);            % relative displacement
xd = y(:,2);           % relative velocity
% Relative acceleration from equation of motion
xdd = -2*zeta*w*xd - w*w*x - ag;
% Absolute acceleration
abs_acc = xdd + ag;
Sa = max(abs(abs_acc));
end
