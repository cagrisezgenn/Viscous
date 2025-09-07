function [geom, orf, hyd, therm] = init_damper_params(geom, orf, hyd, therm)
%INIT_DAMPER_PARAMS Populate damper parameter structs with defaults and derived fields.
%
%   [geom, orf, hyd, therm] = INIT_DAMPER_PARAMS(geom, orf, hyd, therm)
%   fills missing fields in the input structs, computes derived quantities
%   such as piston/orifice areas (Ap, Ao), hydraulic inductance Lh,
%   oil volumes and thermal capacities.
%
%   The number of stories can be provided via hyd.nStories. If absent, a
%   value of 1 is assumed.  hyd.n_parallel controls parallel dampers per
%   story and defaults to 2 when unspecified.
%
%   Inputs:
%       geom  - geometry struct with raw fields (Dp, Lgap, d_o, Lori, etc.)
%       orf   - orifice struct (may be incomplete)
%       hyd   - hydraulic struct (may be incomplete)
%       therm - thermal/temperature struct (may be incomplete)
%
%   Outputs are the same structs augmented with derived fields:
%       geom.Ap, geom.Ap_eff
%       orf.Ao, orf.Ao_eff
%       hyd.Lh, hyd.V0, hyd.V_oil_per, hyd.n_parallel
%       therm.C_oil, therm.C_steel
%
%   See also GETFIELD_DEFAULT.

if nargin < 1 || ~isstruct(geom),  geom  = struct(); end
if nargin < 2 || ~isstruct(orf),   orf   = struct(); end
if nargin < 3 || ~isstruct(hyd),   hyd   = struct(); end
if nargin < 4 || ~isstruct(therm), therm = struct(); end

% --- Orifice defaults ---
orf.n_orf  = getfield_default(orf, 'n_orf',  6);
orf.Cd0    = getfield_default(orf, 'Cd0',    0.61);
orf.CdInf  = getfield_default(orf, 'CdInf',  0.8);
orf.Rec    = getfield_default(orf, 'Rec',    3000);
orf.p_exp  = getfield_default(orf, 'p_exp',  1.1);
orf.p_amb  = getfield_default(orf, 'p_amb',  1.0e5);
orf.cav_sf = getfield_default(orf, 'cav_sf', 0.9);

% --- Thermal defaults ---
therm.antoine_A = getfield_default(therm,'antoine_A',5.0);
therm.antoine_B = getfield_default(therm,'antoine_B',1700);
therm.antoine_C = getfield_default(therm,'antoine_C',-80);
therm.T0_C      = getfield_default(therm,'T0_C',25);
therm.Ts0_C     = getfield_default(therm,'Ts0_C',25);
therm.T_env_C   = getfield_default(therm,'T_env_C',25);
therm.hA_o_env  = getfield_default(therm,'hA_o_env',800);
therm.hA_s_env  = getfield_default(therm,'hA_s_env',600);
therm.hA_os     = getfield_default(therm,'hA_os',   800);
therm.resFactor = getfield_default(therm,'resFactor',22);
therm.cp_oil    = getfield_default(therm,'cp_oil',   1800);
therm.cp_steel  = getfield_default(therm,'cp_steel', 500);
therm.rho_ref   = getfield_default(therm,'rho_ref',  850);
therm.T_ref_C   = getfield_default(therm,'T_ref_C',   25);
therm.alpha_rho = getfield_default(therm,'alpha_rho',7e-4);
therm.beta0     = getfield_default(therm,'beta0',   1.6e9);
therm.b_beta    = getfield_default(therm,'b_beta', -0.0035);
therm.mu_ref    = getfield_default(therm,'mu_ref',   0.9);
therm.b_mu      = getfield_default(therm,'b_mu',   -0.011);

% --- Hydraulic defaults ---
hyd.n_parallel = getfield_default(hyd,'n_parallel',2);
hyd.K_leak     = getfield_default(hyd,'K_leak',     1e-8);
hyd.Vmin_fac   = getfield_default(hyd,'Vmin_fac',0.97);
if ~isfield(hyd,'nStories') || ~isfinite(hyd.nStories) || hyd.nStories<=0
    hyd.nStories = 1;
end

% --- Derived geometric quantities ---
geom.Ap   = pi*geom.Dp^2/4;
orf.Ao    = orf.n_orf * (pi*geom.d_o^2/4);
hyd.Lh    = therm.rho_ref * geom.Lori / ( (orf.n_orf * (pi*geom.d_o^2/4))^2 );

nd = max(1, hyd.n_parallel);
geom.Ap_eff = nd * geom.Ap;
orf.Ao_eff  = nd * orf.Ao;
hyd.n_parallel = nd; % enforce integer

% --- Volumes and heat capacities ---
hyd.V0        = 1.2 * 0.5 * (geom.Ap * (2*geom.Lgap));
hyd.V_oil_per = therm.resFactor * (geom.Ap * (2*geom.Lgap));
steel_to_oil_mass_ratio = 1.5;
nDtot = hyd.nStories * nd;

m_oil_tot   = nDtot * (therm.rho_ref * hyd.V_oil_per);
m_steel_tot = steel_to_oil_mass_ratio * m_oil_tot;
therm.C_oil   = m_oil_tot   * therm.cp_oil;
therm.C_steel = m_steel_tot * therm.cp_steel;

end
