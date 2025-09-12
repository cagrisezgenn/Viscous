function w = pf_weight(t, cfg)
%PF_WEIGHT Time-dependent pressure-force ramp factor.
%   w = PF_WEIGHT(t, cfg) returns the ramp weight for pressure-force
%   activation. It supports scalar or vector inputs t. The weight ramps
%   from 0 to 1 starting at cfg.PF.t_on with time constant cfg.PF.tau.
%
%   The result is multiplied by cfg.on.pressure_force so that the weight
%   is zero when the pressure-force mechanism is disabled. The softness
%   parameter cfg.PF.k controls the smoothness of the ramp onset and the
%   lower bound on cfg.PF.tau via a softplus function.

    k = cfg.PF.k;
    tau_floor = 1e-6;
    dt  = softplus((t - cfg.PF.t_on) ./ k) * k;
    tau = softplus((cfg.PF.tau - tau_floor) ./ k) * k + tau_floor;
    w = cfg.on.pressure_force * (1 - exp(-dt ./ tau));
end

function y = softplus(x)
    y = log1p(exp(-abs(x))) + max(x, 0);
end

