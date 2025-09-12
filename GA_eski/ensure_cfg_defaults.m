function cfg = ensure_cfg_defaults(cfg)
    if nargin==0 || ~isstruct(cfg), cfg = struct(); end
    % top-level toggles
    def_top = {'use_orifice',true; 'use_thermal',true};
    for i=1:size(def_top,1)
        if ~isfield(cfg,def_top{i,1}), cfg.(def_top{i,1}) = def_top{i,2}; end
    end
    % nested on-struct
    if ~isfield(cfg,'on') || ~isstruct(cfg.on), cfg.on = struct(); end
    def_on = {'CdRe',true; 'Rlam',true; 'Rkv',true; 'Qsat',true; 'cavitation',true; ...
              'dP_cap',true; 'hyd_inertia',true; 'leak',true; ...
              'pressure_ode',true; 'pressure_force',true; 'mu_floor',true; ...
              'pf_resistive_only', false};   % <<< YENÄ°: resistive-only anahtarÄ±
    for i=1:size(def_on,1)
        f = def_on{i,1}; if ~isfield(cfg.on,f), cfg.on.(f) = def_on{i,2}; end
    end
    % PF struct
    if ~isfield(cfg,'PF') || ~isstruct(cfg.PF), cfg.PF = struct(); end
    % auto_t_on=true -> t_on = t5 + 0.5 via set_pf_ton_if_auto
    cfg.PF.mode      = getfield_default(cfg.PF,'mode','ramp');
    cfg.PF.auto_t_on = getfield_default(cfg.PF,'auto_t_on', true);
    cfg.PF.t_on      = getfield_default(cfg.PF,'t_on', 0);
    cfg.PF.tau       = getfield_default(cfg.PF,'tau',  2.5);
    cfg.PF.k        = getfield_default(cfg.PF,'k',   0.01);
    cfg.PF.gain      = getfield_default(cfg.PF,'gain', 0.6);
    cfg.PF.gain      = max(0.4, min(0.8, cfg.PF.gain));
    % PF resistive-only clamp options
    % Prefer PF.resistive_only if present; otherwise, mirror legacy cfg.on.pf_resistive_only
    legacy_ro = false; if isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only'), legacy_ro = logical(cfg.on.pf_resistive_only); end
    cfg.PF.resistive_only  = getfield_default(cfg.PF,'resistive_only', legacy_ro);
    % Slope for tanh clamp: s = tanh(resistive_slope * dvel)
    cfg.PF.resistive_slope = getfield_default(cfg.PF,'resistive_slope', 20);
    % Solver options (ODE tolerances and steps)
    if ~isfield(cfg,'solver') || ~isstruct(cfg.solver), cfg.solver = struct(); end
    % Preset controls overall defaults; can be 'explore'|'balanced'|'tight'
    cfg.solver.preset       = getfield_default(cfg.solver, 'preset', 'balanced');
    % Direct overrides (leave empty/NaN to use preset):
    % - RelTol: scalar
    % - AbsTol: scalar (acts as scale factor) or full-length vector (applied as-is)
    % - AbsTolScale: extra multiplier on AbsTol (default 1.0)
    % - MaxStep, InitialStep: seconds
    cfg.solver.RelTol       = getfield_default(cfg.solver, 'RelTol', NaN);
    cfg.solver.AbsTol       = getfield_default(cfg.solver, 'AbsTol', []);
    cfg.solver.AbsTolScale  = getfield_default(cfg.solver, 'AbsTolScale', 1.0);
    cfg.solver.MaxStep      = getfield_default(cfg.solver, 'MaxStep', NaN);
    cfg.solver.InitialStep  = getfield_default(cfg.solver, 'InitialStep', NaN);
    % constraint defaults
    if ~isfield(cfg,'cons') || ~isstruct(cfg.cons), cfg.cons = struct(); end
    if ~isfield(cfg.cons,'spring') || ~isstruct(cfg.cons.spring), cfg.cons.spring = struct(); end
    if ~isfield(cfg.cons.spring,'use_fixed_length'), cfg.cons.spring.use_fixed_length = false; end
    if ~isfield(cfg.cons.spring,'L_free_fixed'), cfg.cons.spring.L_free_fixed = NaN; end
end

