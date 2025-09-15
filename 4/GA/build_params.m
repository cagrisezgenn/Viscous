function params = build_params(params)
%BUILD_PARAMS Compute derived damper and hydraulic fields once.
%   PARAMS = BUILD_PARAMS(PARAMS) fills in fields such as Ap, k_sd,
%   c_lam0, Qcap_big and c_lam_min based on the fundamental geometry and
%   material properties stored in PARAMS. The input struct is returned with
%   the additional fields populated.

if nargin < 1 || ~isstruct(params)
    params = struct();
end

% Reuse existing utility for core damper quantities
params = Utils.recompute_damper_params(params);

% Large orifice flow cap (per damper, adjusted for parallels)
    if isfield(params,'orf') && isfield(params.orf,'CdInf') && ...
            isfield(params,'Ao') && isfield(params,'rho')
    params.Qcap_big = max(params.orf.CdInf * params.Ao, 1e-9) * ...
        sqrt(2*1.0e9 / params.rho);
end

% Minimum laminar damping based on c_lam0
if isfield(params,'c_lam_min_abs') && isfield(params,'c_lam_min_frac') && ...
        isfield(params,'c_lam0')
    params.c_lam_min = max(params.c_lam_min_abs, ...
        params.c_lam_min_frac * params.c_lam0);
end

end
