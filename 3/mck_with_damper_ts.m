function [x,a_rel,ts,diag] = mck_with_damper_ts(t,ag,M,C,K,k_sd,c_lam0,use_orifice,orf,rho,Ap,Ao,Qcap,mu_ref,use_thermal,thermal,T0_C,T_ref_C,b_mu,c_lam_min,c_lam_cap,Lgap,cp_oil,cp_steel,steel_to_oil_mass_ratio,toggle_gain,story_mask,n_dampers_per_story,resFactor,cfg,F_story_target)
%MCK_WITH_DAMPER_TS Wrapper around MCK_WITH_DAMPER returning time-series.
%
%   [X,A_REL,TS,DIAG] = MCK_WITH_DAMPER_TS(T,AG,M,C,K, ...) calls the existing
%   MCK_WITH_DAMPER function to compute the structural response and then
%   assembles a time-series structure TS containing various per–time-step
%   diagnostics such as power and energy histories. The original diagnostic
%   structure DIAG is returned unchanged.

if nargin < 31 || isempty(F_story_target)
    F_story_target = [];
end

[x,a_rel,diag] = mck_with_damper(t,ag,M,C,K,k_sd,c_lam0,use_orifice,orf,rho,Ap,Ao,Qcap,mu_ref,...
    use_thermal,thermal,T0_C,T_ref_C,b_mu,c_lam_min,c_lam_cap,Lgap,cp_oil,cp_steel,steel_to_oil_mass_ratio,
    toggle_gain,story_mask,n_dampers_per_story,resFactor,cfg,F_story_target);

% Story vectors needed for power calculations
nStories = size(diag.drift,2);
Rvec = toggle_gain(:); if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
mask = story_mask(:); if numel(mask)==1, mask = mask*ones(nStories,1); end
ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
multi = (mask .* ndps).';
Rvec = Rvec.';

% Basic time series fields from diagnostics
ts = struct();
ts.t = t;
ts.drift = diag.drift;
ts.dvel = diag.dvel;
ts.story_force = diag.story_force;
ts.PF = diag.PF;
ts.Q = diag.Q;
ts.dP_orf = diag.dP_orf;
ts.Qcap_ratio = abs(diag.Q) ./ Qcap;
ts.cav_mask = diag.dP_orf < 0;

% Power components
P_visc_per = diag.c_lam .* (diag.dvel.^2);
ts.P_visc = sum(P_visc_per .* multi, 2);
P_orf_per = diag.dP_orf .* diag.Q;
ts.P_orf = sum(P_orf_per .* multi, 2);
if isfield(diag,'P_sum')
    ts.P_sum = diag.P_sum;
elseif isfield(ts,'P_orf') && isfield(ts,'P_visc')
    ts.P_sum = ts.P_orf + ts.P_visc;
else
    ts.P_sum = [];
end

% Energy accumulations
P_struct = sum(diag.story_force .* diag.dvel, 2);
ts.E_orf = cumtrapz(t, ts.P_orf);
ts.E_struct = cumtrapz(t, P_struct);
% mechanical energy integral not exported/consumed

% DIAG is returned unchanged
end
