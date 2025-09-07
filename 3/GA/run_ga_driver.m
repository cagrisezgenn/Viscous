function [X,F,gaout] = run_ga_driver(scaledOrSnap, params, optsEval, optsGA)
% === SAFE PARPOOL OPEN (cleanup + cap threads) ===
try
    parpool_hard_reset(16);
catch
end
%RUN_GA_DRIVER Hybrid GA driver: accepts snapshot path or in-memory structs.
%   [X,F,GAOUT] = RUN_GA_DRIVER(SCALED_OR_PATH, PARAMS, OPTSEVAL, OPTSGA)
%   Two modes:
%     A) Snapshot: run_ga_driver('out/<ts>/snapshot.mat', params?, ...)
%        Loads 'scaled' (and params if not provided).
%     B) In-memory: run_ga_driver(scaled, params, ...)
%   Fitness evaluations have no IO; no re-scaling is performed.

    % ---------- Zero-arg convenience: pull from base workspace ----------
    % Avoid referencing missing input args when nargin==0/1 by using locals.
    if nargin >= 1
        scaledOrSnap_local = scaledOrSnap;
    else
        scaledOrSnap_local = [];
    end
    if nargin >= 2
        params_local = params;
    else
        params_local = [];
    end
    % If called without inputs (e.g., pressing Run in editor), try to fetch
    % 'scaled' and 'params' from the base workspace so that GA can start.
    if isempty(scaledOrSnap_local)
        try
            scaledOrSnap_local = evalin('base','scaled');
        catch
            % leave empty; handled by asserts below
        end
    end
    if isempty(params_local)
        try
            params_local = evalin('base','params');
        catch
            % may be filled from snapshot path below
        end
    end
    if nargin < 3 || isempty(optsEval), optsEval = struct; end
    if nargin < 4 || isempty(optsGA),   optsGA   = struct; end

    % ---------- Resolve scaled/params/meta ----------
    if ischar(scaledOrSnap_local) || isstring(scaledOrSnap_local)
        S = load(char(scaledOrSnap_local), 'scaled','params','opts', ...
                 'IM_mode','band_fac','s_bounds','mu_factors','mu_weights','thr');
        assert(isfield(S,'scaled') && ~isempty(S.scaled), ...
            'run_ga_driver: snapshot missing ''scaled''.');
        scaled = S.scaled;
        if isempty(params_local)
            assert(isfield(S,'params') && ~isempty(S.params), ...
                'run_ga_driver: params not supplied and snapshot missing params.');
            params = S.params;
        else
            params = params_local;
        end
        meta = struct();
        meta.IM_mode    = Utils.getfield_default(S,'IM_mode','');
        meta.band_fac   = Utils.getfield_default(S,'band_fac',[]);
        meta.s_bounds   = Utils.getfield_default(S,'s_bounds',[]);
        meta.mu_factors = Utils.getfield_default(S,'mu_factors',[0.75 1.00 1.25]);
        meta.mu_weights = Utils.getfield_default(S,'mu_weights',[0.2 0.6 0.2]);
        thr_default = struct('dP95_max',50e6,'Qcap95_max',0.5,'cav_pct_max',0,'T_end_max',75,'mu_end_min',0.5);
        meta.thr       = Utils.getfield_default(S,'thr', thr_default);
    else
        scaled = scaledOrSnap_local;
        params = params_local;
        thr_default = struct('dP95_max',50e6,'Qcap95_max',0.5,'cav_pct_max',0,'T_end_max',75,'mu_end_min',0.5);
        meta = struct('IM_mode','', 'band_fac',[], 's_bounds',[], ...
                      'mu_factors',[0.75 1.00 1.25], 'mu_weights',[0.2 0.6 0.2], ...
                      'thr', thr_default);
        % Auto-prepare workspace if needed (no inputs provided)
        if (isempty(scaled) || isempty(params))
            try
                % Ensure paths
                try, setup; catch, end
                % Load base parameters and compute T1
                parametreler;
                % Freeze dataset (band scaling + caps + TRIM); no export
                try
                    [~, scaled] = load_ground_motions(T1);
                catch
                    % fall back to raw if needed
                    [scaled_raw, ~] = load_ground_motions(); %#ok<ASGLU>
                    error('run_ga_driver:auto_prep','Failed to produce scaled set.');
                end
                % Build params struct (mirror damperlinon)
                params = struct('M',M,'C0',C0,'K',K,'k_sd',k_sd,'c_lam0',c_lam0, ...
                    'orf',orf,'rho',rho,'Ap',Ap,'A_o',A_o,'Qcap_big',Qcap_big,'mu_ref',mu_ref, ...
                    'thermal',thermal,'T0_C',T0_C,'T_ref_C',T_ref_C,'b_mu',b_mu, ...
                    'c_lam_min',c_lam_min,'c_lam_cap',c_lam_cap,'Lgap',Lgap, ...
                    'cp_oil',cp_oil,'cp_steel',cp_steel,'steel_to_oil_mass_ratio',steel_to_oil_mass_ratio, ...
                    'toggle_gain',toggle_gain,'story_mask',story_mask,'n_dampers_per_story',n_dampers_per_story, ...
                    'resFactor',resFactor,'cfg',cfg,'story_height',story_height);
                % Export to base workspace for user convenience
                try
                    assignin('base','scaled',scaled);
                    assignin('base','params',params);
                    assignin('base','T1',T1);
                catch
                end
            catch
                % leave to asserts below
            end
        end
    end
    assert(~isempty(scaled),'run_ga_driver: scaled set is empty. Define ''scaled'' in workspace or pass a snapshot path.');
    assert(~isempty(params),'run_ga_driver: params is empty. Define ''params'' in workspace or include it in snapshot.');

    % ---------- Defaults for evaluation (IO-free) ----------
    if nargin < 3 || isempty(optsEval), optsEval = struct; end
    optsEval.do_export     = false;
    optsEval.quiet         = true;
    optsEval.thermal_reset = 'each';
    if ~isfield(optsEval,'mu_factors'), optsEval.mu_factors = meta.mu_factors; end
    if ~isfield(optsEval,'mu_weights'), optsEval.mu_weights = meta.mu_weights; end
    if ~isfield(optsEval,'thr'),        optsEval.thr        = meta.thr; end


    % ---------- GA options and objective ----------
    if nargin < 4 || isempty(optsGA), optsGA = struct; end
    rng(42);

    % decision vector: [d_o_mm, n_orf, g_lo, g_mid, g_hi, PF_tau, PF_gain]
    lb = [2.80, 5, 3.60, 3.80, 1.50, 0.95, 0.78];
    ub = [3.60, 6, 4.00, 4.00, 3.60, 1.10, 0.90];
    IntCon = 2;  % n_orf only

    % Provide dataset signature for memo keys
    try, dsig = sum([scaled.IM]) + sum([scaled.PGA]); catch, dsig = 0; end
    optsEval.dsig = dsig;
    obj = @(x) eval_design_fast(x, scaled, params, optsEval); % quantize/clamp inside

    options = optimoptions('gamultiobj', ...
       'PopulationSize',    Utils.getfield_default(optsGA,'PopulationSize',16), ...
       'MaxGenerations',    Utils.getfield_default(optsGA,'MaxGenerations',8), ...
       'CrossoverFraction', Utils.getfield_default(optsGA,'CrossoverFraction',0.85), ...
       'MutationFcn',       Utils.getfield_default(optsGA,'MutationFcn',{@mutationadaptfeasible}), ...
       'ParetoFraction',    Utils.getfield_default(optsGA,'ParetoFraction',0.70), ...
       'StallGenLimit',     Utils.getfield_default(optsGA,'StallGenLimit',10), ...
       'DistanceMeasureFcn','distancecrowding', ...
       'UseParallel',       Utils.getfield_default(optsGA,'UseParallel',true), ...
       'Display','iter','PlotFcn',[], 'FunctionTolerance',1e-3);

    % Provide grid-aligned initial population if none supplied (sparser grid + seeds)
    try
        if ~isfield(optsGA,'InitialPopulationMatrix') || isempty(optsGA.InitialPopulationMatrix)
            step_vec = [0.1 NaN 0.05 0.05 0.05 0.05 0.02];
            P0 = Utils.initial_pop_grid(lb, ub, options.PopulationSize, step_vec);
            % Seeds (feasible within new bounds)
            seed = [ 2.80 5 3.60 3.80 1.50 0.95 0.78;
                     3.20 5 3.80 3.90 2.50 1.00 0.82;
                     3.60 6 4.00 4.00 3.60 1.10 0.90 ];
            ns = min(size(seed,1), size(P0,1));
            P0(1:ns,:) = seed(1:ns,:);
            % Try to read latest ga_front.csv to seed with prior Pareto
            try
                dd = dir(fullfile('out','ga_*'));
                if ~isempty(dd)
                    [~,ix] = max([dd.datenum]);
                    latest = fullfile(dd(ix).folder, dd(ix).name, 'ga_front.csv');
                    if exist(latest,'file')
                        Tprev = readtable(latest);
                        cols = {'d_o_mm','n_orf','g_lo','g_mid','g_hi','PF_tau','PF_gain'};
                        if all(ismember(cols, Tprev.Properties.VariableNames))
                            Xprev = Tprev{:, cols};
                            % clamp & quantize to current grid
                            for irow = 1:size(Xprev,1)
                                Xprev(irow,:) = quant_clamp_x(Xprev(irow,:));
                            end
                            % mix: previous + P0; keep unique and fit to pop size
                            Pmix = unique([Xprev; P0], 'rows', 'stable');
                            if size(Pmix,1) >= options.PopulationSize
                                P0 = Pmix(1:options.PopulationSize,:);
                            else
                                need = options.PopulationSize - size(Pmix,1);
                                P0 = [Pmix; P0(1:need,:)];
                            end
                        end
                    end
                end
            catch
            end
            options = optimoptions(options,'InitialPopulationMatrix', P0);
        end
    catch
    end

    [X,F,exitflag,output,population,scores] = gamultiobj(obj, numel(lb), [],[],[],[], lb, ub, [], IntCon, options);
    gaout = struct('exitflag',exitflag,'output',output);

    % ---------- Post-run lightweight packaging (no simulations) ----------
    outdir = fullfile('out', ['ga_' datestr(now,'yyyymmdd_HHMMSS_FFF')]);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    opts_ga = options; date_str = datestr(now);
    save(fullfile(outdir,'ga_front.mat'),'X','F','opts_ga','meta','date_str','-v7.3');

    % === Re-evaluate Pareto designs to collect metrics per-row ===
    nF = size(X,1);
    x10pk  = zeros(nF,1);  a10pk = zeros(nF,1);
    dP95   = zeros(nF,1);  Qcap95 = zeros(nF,1); cavW = zeros(nF,1);
    Tend   = zeros(nF,1);  muend  = zeros(nF,1);
    PFp95  = zeros(nF,1);  Qq50 = zeros(nF,1);  Qq95 = zeros(nF,1);
    dPq50  = zeros(nF,1);  dPq95w = zeros(nF,1); Toil = zeros(nF,1); Tsteel = zeros(nF,1);
    Etot   = zeros(nF,1);  Eor = zeros(nF,1);   Estr  = zeros(nF,1); Eratio = zeros(nF,1); Pmech = zeros(nF,1);

    % penalty parts (same as eval)
    pen     = zeros(nF,1);
    pen_dP  = zeros(nF,1); pen_Qcap = zeros(nF,1); pen_cav = zeros(nF,1); pen_T = zeros(nF,1); pen_mu = zeros(nF,1);
    lambda  = Utils.getfield_default(optsEval,'penalty_scale',10);
    pwr     = Utils.getfield_default(optsEval,'penalty_power',1.0);
    W       = Utils.getfield_default(optsEval,'penalty_weights', struct('dP',1,'Qcap',1,'cav',2,'T',0.5,'mu',0.5));
    rel = @(v,lim) max(0,(v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0,(lim - v)./max(lim,eps)).^pwr;

    Opost = struct('do_export',false,'quiet',true,'thermal_reset','each','order','natural', ...
                   'use_orifice', true, 'use_thermal', true, ...
                   'mu_factors', meta.mu_factors, 'mu_weights', meta.mu_weights, 'thr', meta.thr);

    parfor i = 1:nF
        Xi = quant_clamp_x(X(i,:));
        Pi = decode_params_from_x(params, Xi);
        Si = run_batch_windowed(scaled, Pi, Opost);

        % Extract table columns safely
        try, v_x10 = Si.table.x10_max_D_worst; catch, v_x10 = 0; end
        try, v_a10 = Si.table.a10abs_max_D_worst; catch, v_a10 = 0; end

        try, v_dP95 = Si.table.dP95_worst; catch, v_dP95 = 0; end
        try, v_Qcap = Si.table.Qcap95_worst; catch, v_Qcap = 0; end
        try, v_cav  = Si.table.cav_pct_worst; catch, v_cav = 0; end
        try, v_Tend = Si.table.T_end_worst; catch, v_Tend = 0; end
        try, v_mu   = Si.table.mu_end_worst; catch, v_mu = 1; end

        try, v_PFp95 = Si.table.PF_p95_worst; catch, v_PFp95 = 0; end
        try, v_Qq50  = Si.table.Q_q50_worst; catch, v_Qq50 = 0; end
        try, v_Qq95  = Si.table.Q_q95_worst; catch, v_Qq95 = 0; end
        try, v_dPq50 = Si.table.dP_orf_q50_worst; catch, v_dPq50 = 0; end
        try, v_Toil  = Si.table.T_oil_end_worst; catch, v_Toil = []; end
        try, v_Tsteel= Si.table.T_steel_end_worst; catch, v_Tsteel = []; end

        try, v_Eor = Si.table.E_orifice_sum; catch, v_Eor = 0; end
        try, v_Est = Si.table.E_struct_sum;  catch, v_Est = 0; end
        try, v_Pm  = Si.table.P_mech_sum;    catch, v_Pm = 0; end

        % Aggregate across records (dataset) to scalars per design
        x10pk(i)  = max(v_x10(:));
        a10pk(i)  = max(v_a10(:));

        dP95(i)   = max(v_dP95(:));
        Qcap95(i) = max(v_Qcap(:));
        cavW(i)   = max(v_cav(:));
        Tend(i)   = max(v_Tend(:));
        if isempty(v_mu), v_mu = 1; end
        muend(i)  = min(v_mu(:));

        PFp95(i)  = max(v_PFp95(:));
        Qq50(i)   = max(v_Qq50(:));
        Qq95(i)   = max(v_Qq95(:));
        dPq50(i)  = max(v_dPq50(:));
        dPq95w(i) = max(v_dP95(:));
        if isempty(v_Toil), v_Toil = v_Tend; end
        Toil(i)   = max(v_Toil(:));
        if isempty(v_Tsteel), v_Tsteel = 0; end
        Tsteel(i) = max(v_Tsteel(:));

        Eor(i)    = sum(v_Eor(:));
        Estr(i)   = sum(v_Est(:));
        Etot(i)   = Eor(i) + Estr(i);
        Eratio(i) = (Estr(i)>0) * (Eor(i)/max(Estr(i),eps));
        Pmech(i)  = sum(v_Pm(:));

        % Penalty parts (mirror eval_design_fast): mean over records
        vv_dP   = v_dP95(:); if isempty(vv_dP), vv_dP = 0; end
        vv_Qcap = v_Qcap(:); if isempty(vv_Qcap), vv_Qcap = 0; end
        vv_cav  = v_cav(:);  if isempty(vv_cav), vv_cav = 0; end
        vv_Tend = v_Tend(:); if isempty(vv_Tend), vv_Tend = 0; end
        vv_mu   = v_mu(:);   if isempty(vv_mu), vv_mu = 1; end

        pen_dP(i)   = mean(rel(vv_dP,   Opost.thr.dP95_max));
        pen_Qcap(i) = mean(rel(vv_Qcap, Opost.thr.Qcap95_max));
        if Opost.thr.cav_pct_max<=0
            pen_cav(i) = mean(max(0,vv_cav).^pwr);
        else
            pen_cav(i) = mean(rel(vv_cav, Opost.thr.cav_pct_max));
        end
        pen_T(i)    = mean(rel(vv_Tend, Opost.thr.T_end_max));
        pen_mu(i)   = mean(rev(vv_mu,   Opost.thr.mu_end_min));
        pen(i)      = lambda*(W.dP*pen_dP(i)+W.Qcap*pen_Qcap(i)+W.cav*pen_cav(i)+W.T*pen_T(i)+W.mu*pen_mu(i));
    end

    % Build T using per-row arrays
    T = array2table([X F pen pen_dP pen_Qcap pen_cav pen_T pen_mu], 'VariableNames', ...
       {'d_o_mm','n_orf','g_lo','g_mid','g_hi','PF_tau','PF_gain','f1','f2', ...
        'pen','pen_dP','pen_Qcap','pen_cav','pen_T','pen_mu'});

    T.x10_max_damperli    = x10pk;
    T.a10abs_max_damperli = a10pk;
    % extra diagnostics wanted
    T.dP95_worst          = dP95;     T.Qcap95_worst      = Qcap95;   T.cav_pct_worst = cavW;
    T.T_end_worst         = Tend;     T.mu_end_worst      = muend;    T.PF_p95_worst  = PFp95;
    T.Q_q50_worst         = Qq50;     T.Q_q95_worst       = Qq95;     T.dP_orf_q50_worst = dPq50;
    T.dP_orf_q95_worst    = dPq95w;   T.T_oil_end_worst   = Toil;     T.T_steel_end_worst = Tsteel;
    T.energy_tot_sum      = Etot;     T.E_orifice_sum     = Eor;      T.E_struct_sum  = Estr;
    T.E_ratio             = Eratio;   T.P_mech_sum        = Pmech;

    % === BASELINE (pre-GA) ROW: params başlangıcıyla tek koşu, ilk satır ===
    try
        % Build X0 from params (with sensible fallbacks)
        X0 = nan(1,7);
        % d_o_mm and n_orf
        try
            if isfield(params,'orf') && isfield(params.orf,'d_o') && ~isempty(params.orf.d_o)
                X0(1) = 1e3 * params.orf.d_o;
            else
                X0(1) = 3.00;
            end
        catch, X0(1) = 3.00; end
        try
            if isfield(params,'n_orf') && ~isempty(params.n_orf)
                X0(2) = params.n_orf;
            else
                X0(2) = 5;
            end
        catch, X0(2) = 5; end
        % g_lo, g_mid, g_hi from toggle_gain if available
        try
            if isfield(params,'toggle_gain') && ~isempty(params.toggle_gain)
                tg = params.toggle_gain(:);
                nStories = numel(tg);
                loN = min(3, nStories); hiN = min(2, nStories);
                midIdx = (loN+1) : max(nStories-hiN, loN+1);
                if isempty(midIdx), midIdx = round(nStories/2); end
                X0(3) = mean(tg(1:loN));
                X0(4) = mean(tg(midIdx));
                X0(5) = mean(tg(end-hiN+1:end));
            else
                X0(3:5) = [3.90 3.95 1.50];
            end
        catch, X0(3:5) = [3.90 3.95 1.50]; end
        % PF parameters
        try
            if isfield(params,'cfg') && isfield(params.cfg,'PF')
                X0(6) = Utils.getfield_default(params.cfg.PF,'tau',0.97);
                X0(7) = Utils.getfield_default(params.cfg.PF,'gain',0.90);
            else
                X0(6:7) = [0.97 0.90];
            end
        catch, X0(6:7) = [0.97 0.90]; end

        % Simulate baseline with same post-eval options
        X0 = quant_clamp_x(X0);
        P0    = decode_params_from_x(params, X0);
        S0    = run_batch_windowed(scaled, P0, Opost);
        T0bl  = S0.table;
        % Objectives
        f0 = [NaN NaN];
        try
            [f0, ~] = eval_design_fast(X0, scaled, params, Opost);
        catch
        end
        % Penalty parts (same formula)
        dP95_0   = max(T0bl.dP95_worst);
        Qcap95_0 = max(T0bl.Qcap95_worst);
        cavW_0   = max(T0bl.cav_pct_worst);
        Tend_0   = max(T0bl.T_end_worst);
        muend_0  = min(T0bl.mu_end_worst);
        pen_dP_0   = rel(dP95_0,  Opost.thr.dP95_max);
        pen_Qcap_0 = rel(Qcap95_0,Opost.thr.Qcap95_max);
        if Opost.thr.cav_pct_max<=0, pen_cav_0 = max(0,cavW_0).^pwr; else, pen_cav_0 = rel(cavW_0,Opost.thr.cav_pct_max); end
        pen_T_0    = rel(Tend_0,  Opost.thr.T_end_max);
        pen_mu_0   = rev(muend_0, Opost.thr.mu_end_min);
        pen_0      = lambda*(W.dP*pen_dP_0 + W.Qcap*pen_Qcap_0 + W.cav*pen_cav_0 + W.T*pen_T_0 + W.mu*pen_mu_0);

        % Build baseline row T0 with same columns as T
        vn = T.Properties.VariableNames;
        T0 = cell2table(cell(1,numel(vn)), 'VariableNames', vn);
        % decision variables
        T0.d_o_mm = X0(1); T0.n_orf = X0(2); T0.g_lo = X0(3); T0.g_mid = X0(4); T0.g_hi = X0(5);
        T0.PF_tau = X0(6); T0.PF_gain = X0(7);
        % objectives
        T0.f1 = f0(1); T0.f2 = f0(2);
        % penalties
        T0.pen = pen_0; T0.pen_dP = pen_dP_0; T0.pen_Qcap = pen_Qcap_0; T0.pen_cav = pen_cav_0; T0.pen_T = pen_T_0; T0.pen_mu = pen_mu_0;
        % damper peak
        try
            if ismember('x10_max_damperli', vn)
                T0.x10_max_damperli = max(T0bl.x10_max_D_worst);
            end
            if ismember('a10abs_max_damperli', vn)
                T0.a10abs_max_damperli = max(T0bl.a10abs_max_D_worst);
            end
        catch
        end
        % other diagnostics if present
        try
            if ismember('dP95_worst', vn),        T0.dP95_worst        = max(T0bl.dP95_worst); end
            if ismember('Qcap95_worst', vn),      T0.Qcap95_worst      = max(T0bl.Qcap95_worst); end
            if ismember('cav_pct_worst', vn),     T0.cav_pct_worst     = max(T0bl.cav_pct_worst); end
            if ismember('T_end_worst', vn),       T0.T_end_worst       = max(T0bl.T_end_worst); end
            if ismember('mu_end_worst', vn),      T0.mu_end_worst      = min(T0bl.mu_end_worst); end
            if ismember('PF_p95_worst', vn) && ismember('PF_p95_worst', T0bl.Properties.VariableNames)
                T0.PF_p95_worst = max(T0bl.PF_p95_worst);
            end
            if ismember('Q_q50_worst', vn) && ismember('Q_q50_worst', T0bl.Properties.VariableNames)
                T0.Q_q50_worst = max(T0bl.Q_q50_worst);
            end
            if ismember('Q_q95_worst', vn) && ismember('Q_q95_worst', T0bl.Properties.VariableNames)
                T0.Q_q95_worst = max(T0bl.Q_q95_worst);
            end
            if ismember('dP_orf_q50_worst', vn) && ismember('dP_orf_q50_worst', T0bl.Properties.VariableNames)
                T0.dP_orf_q50_worst = max(T0bl.dP_orf_q50_worst);
            end
            if ismember('dP_orf_q95_worst', vn)
                if ismember('dP_orf_q95_worst', T0bl.Properties.VariableNames)
                    T0.dP_orf_q95_worst = max(T0bl.dP_orf_q95_worst);
                else
                    T0.dP_orf_q95_worst = max(T0bl.dP95_worst);
                end
            end
            if ismember('T_oil_end_worst', vn) && ismember('T_oil_end_worst', T0bl.Properties.VariableNames)
                T0.T_oil_end_worst = max(T0bl.T_oil_end_worst);
            end
            if ismember('T_steel_end_worst', vn) && ismember('T_steel_end_worst', T0bl.Properties.VariableNames)
                T0.T_steel_end_worst = max(T0bl.T_steel_end_worst);
            end
            if ismember('E_orifice_sum', vn) && ismember('E_orifice_sum', T0bl.Properties.VariableNames)
                T0.E_orifice_sum = sum(T0bl.E_orifice_sum);
            end
            if ismember('E_struct_sum', vn) && ismember('E_struct_sum', T0bl.Properties.VariableNames)
                T0.E_struct_sum = sum(T0bl.E_struct_sum);
            end
            if ismember('energy_tot_sum', vn)
                if ismember('energy_tot_sum', T0bl.Properties.VariableNames)
                    T0.energy_tot_sum = sum(T0bl.energy_tot_sum);
                else
                    try
                        T0.energy_tot_sum = T0.E_orifice_sum + T0.E_struct_sum;
                    catch, T0.energy_tot_sum = 0; end
                end
            end
            if ismember('E_ratio', vn)
                try
                    T0.E_ratio = (T0.E_struct_sum>0) * (T0.E_orifice_sum / max(T0.E_struct_sum, eps));
                catch, T0.E_ratio = 0; end
            end
            if ismember('P_mech_sum', vn) && ismember('P_mech_sum', T0bl.Properties.VariableNames)
                T0.P_mech_sum = sum(T0bl.P_mech_sum);
            end
        catch
        end

        % Prepend baseline row
        T = [T0; T];
    catch
    end

    writetable(T, fullfile(outdir,'ga_front.csv'));

    % Knee-point (min distance to ideal after min-max normalization) with full columns
    try
        f1v = T.f1; f2v = T.f2;
        f1n = (f1v - min(f1v)) / max(eps, (max(f1v)-min(f1v)));
        f2n = (f2v - min(f2v)) / max(eps, (max(f2v)-min(f2v)));
        d   = hypot(f1n, f2n);
        [~,kidx] = min(d);
        Tknee = T(kidx,:);
        try
            Tknee_full = [T(1,:); Tknee]; % assume T(1,:) is baseline just prepended
            writetable(Tknee_full, fullfile(outdir,'ga_knee.csv'));
        catch
            writetable(Tknee, fullfile(outdir,'ga_knee.csv'));
        end
    catch
    end

    % Decode top-K designs without running sims
    K = min(10, size(X,1));
    if K > 0
        idx = unique(round(linspace(1,size(X,1),K)));
        params_list = cell(numel(idx),1);
        for i = 1:numel(idx)
            params_list{i} = decode_params_from_x(params, X(idx(i),:)); %#ok<AGROW>
        end
        save(fullfile(outdir,'ga_front.mat'),'params_list','-append');

        % ga_topK.csv with key decoded fields
        top = table();
        for i = 1:numel(idx)
            P = params_list{i};
            row = table(P.orf.d_o*1e3, P.n_orf, ...
                X(idx(i),3),X(idx(i),4),X(idx(i),5),X(idx(i),6),X(idx(i),7), ...
                P.A_o, P.Qcap_big, 'VariableNames', ...
                {'d_o_mm','n_orf','g_lo','g_mid','g_hi','PF_tau','PF_gain','A_o','Qcap_big'});
            top = [top; row]; %#ok<AGROW>
        end
        writetable(top, fullfile(outdir,'ga_topK.csv'));
    end

end

function [f, meta] = eval_design_fast(x, scaled, params0, optsEval)
    % snap to grids
    x = x(:)';
    % quantize
    x(1) = Utils.quantize_step(x(1),0.05);  % d_o_mm
    x(3) = Utils.quantize_step(x(3),0.05);  % g_lo
    x(4) = Utils.quantize_step(x(4),0.05);  % g_mid
    x(5) = Utils.quantize_step(x(5),0.05);  % g_hi
    x(6) = Utils.quantize_step(x(6),0.01);  % PF_tau
    x(7) = Utils.quantize_step(x(7),0.02);  % PF_gain
    x(2) = round(max(x(2),1));              % n_orf integer, >=1

    % clamp to GA bounds to be safe after quantization
    x(1) = min(max(x(1), 2.80), 3.60);
    x(2) = min(max(x(2), 5), 6);
    x(3) = min(max(x(3), 3.60), 4.00);
    x(4) = min(max(x(4), 3.80), 4.00);
    x(5) = min(max(x(5), 1.50), 3.60);
    x(6) = min(max(x(6), 0.95), 1.10);
    x(7) = min(max(x(7), 0.78), 0.90);

    persistent memo;
    if isempty(memo), memo = containers.Map(); end
    % Add dataset salt so memo keys do not collide across different scaled sets
    dsig = 0;
    try, dsig = sum([scaled.IM]) + sum([scaled.PGA]); catch, dsig = 0; end
    key = jsonencode([x, dsig]);
    if isKey(memo, key)
        meta = memo(key); f = meta.f; return;
    end

    P = decode_params_from_x(params0, x);

    O = struct();
    if nargin >= 4 && ~isempty(optsEval), O = optsEval; end
    O.do_export = false;
    O.quiet  = true;
    O.thermal_reset = 'each';
    O.order = 'natural';
    O.use_orifice = true; O.use_thermal = true;
    if ~isfield(O,'mu_factors'), O.mu_factors = [0.75 1.00 1.25]; end
    if ~isfield(O,'mu_weights'), O.mu_weights = [0.2 0.6 0.2]; end

    % safe evaluation (no IO during GA)
    try
        S = run_batch_windowed(scaled, P, O);
    catch ME
        f = [1e6, 1e6];
        meta = struct('x',x,'f',f,'error','eval_failed', ...
                      'message',ME.message,'identifier',ME.identifier);
        % --- penalties: force numeric (no NaNs) on failure
        meta.pen      = 0;
        meta.pen_dP   = 0;
        meta.pen_Qcap = 0;
        meta.pen_cav  = 0;
        meta.pen_T    = 0;
        meta.pen_mu   = 0;
        try, memo_store('set', jsonencode([x, dsig]), meta); catch, end
        return;
    end

    % --- HARD FILTER (erken eleme) ---
    dP95v = S.table.dP95_worst;
    qcapv = S.table.Qcap95_worst;
    cavv  = S.table.cav_pct_worst;
    if any(dP95v > 1e9) || any(qcapv > 0.90) || any(cavv > 0.01)
        f = [1e6, 1e6];
        meta = struct('x',x,'f',f,'hard_kill',true);
        % --- penalties: ensure numeric zeros for hard-kill
        meta.pen      = 0;
        meta.pen_dP   = 0;
        meta.pen_Qcap = 0;
        meta.pen_cav  = 0;
        meta.pen_T    = 0;
        meta.pen_mu   = 0;
        try, memo_store('set', jsonencode([x, dsig]), meta); catch, end
        return;
    end

    f1 = mean(S.table.PFA_w);
    f2 = mean(S.table.IDR_w);
    thr = O.thr;
    dP95v   = S.table.dP95_worst;    
    qcapv   = S.table.Qcap95_worst;
    cavv    = S.table.cav_pct_worst;
    Tendv   = S.table.T_end_worst;
    muendv  = S.table.mu_end_worst;
    % --- PENALTY (soft, multiplicative) ---
    lambda  = Utils.getfield_default(optsEval,'penalty_scale',10);
    pwr     = Utils.getfield_default(optsEval,'penalty_power',1.0);
    rel = @(v,lim) max(0, (v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0, (lim - v)./max(lim,eps)).^pwr;
    W = Utils.getfield_default(optsEval,'penalty_weights', struct('dP',1,'Qcap',1,'cav',2,'T',0.5,'mu',0.5));
    pen_dP   = isfield(thr,'dP95_max')   * mean(rel(dP95v,  thr.dP95_max));
    pen_Qcap = isfield(thr,'Qcap95_max') * mean(rel(qcapv,  thr.Qcap95_max));
    if isfield(thr,'cav_pct_max')
        if thr.cav_pct_max<=0, pen_cav = mean(max(0,cavv).^pwr); else, pen_cav = mean(rel(cavv, thr.cav_pct_max)); end
    else
        pen_cav = 0;
    end
    pen_T    = isfield(thr,'T_end_max')  * mean(rel(Tendv,  thr.T_end_max));
    pen_mu   = isfield(thr,'mu_end_min') * mean(rev(muendv, thr.mu_end_min));
    pen = W.dP*pen_dP + W.Qcap*pen_Qcap + W.cav*pen_cav + W.T*pen_T + W.mu*pen_mu;
    % multiplicative penalty keeps scale of objectives
    f = [f1, f2] .* (1 + lambda*pen);
    % Build meta and enforce numeric penalties
    meta = struct('x',x,'f',f,'PFA_w_mean',f1,'IDR_w_mean',f2);
    Penalty   = lambda*pen;
    pen_parts = struct('dP',pen_dP,'Qcap',pen_Qcap,'cav',pen_cav,'T',pen_T,'mu',pen_mu);
    % --- penalties: force numeric (no NaNs)
    if ~exist('Penalty','var') || ~isfinite(Penalty), Penalty = 0; end
    if ~exist('pen_parts','var') || ~isstruct(pen_parts), pen_parts = struct(); end
    pf = {'dP','Qcap','cav','T','mu'};
    for ii=1:numel(pf)
        fn = pf{ii};
        if ~isfield(pen_parts,fn) || ~isfinite(pen_parts.(fn)), pen_parts.(fn) = 0; end
    end
    meta.pen      = Penalty;
    meta.pen_dP   = pen_parts.dP;
    meta.pen_Qcap = pen_parts.Qcap;
    meta.pen_cav  = pen_parts.cav;
    meta.pen_T    = pen_parts.T;
    meta.pen_mu   = pen_parts.mu;
    meta.pen_parts = pen_parts;
    % --- append damperli peak metrics (zeros if unavailable)
    x10pk = 0; a10pk = 0;
    try
        if isfield(S,'table') && istable(S.table)
            if ismember('x10_max_D_worst', S.table.Properties.VariableNames)
                try
                    x10pk = max(S.table.x10_max_D_worst);
                catch
                    x10pk = S.table.x10_max_D_worst;
                    if numel(x10pk)>1, x10pk = max(x10pk(:)); end
                end
            end
            if ismember('a10abs_max_D_worst', S.table.Properties.VariableNames)
                try
                    a10pk = max(S.table.a10abs_max_D_worst);
                catch
                    a10pk = S.table.a10abs_max_D_worst;
                    if numel(a10pk)>1, a10pk = max(a10pk(:)); end
                end
            end
        end
    catch
    end
    meta.x10_max_damperli    = x10pk;
    meta.a10abs_max_damperli = a10pk;
    % === Append penalty drivers & diagnostics from S.table (if available) ===
    try
        if isfield(S,'table') && istable(S.table)
            candCols = { ...
              'dP95_worst', 'Qcap95_worst', 'cav_pct_worst', 'T_end_worst', 'mu_end_worst', ...
              'PF_p95_worst', ...
              'Q_q50_worst','Q_q95_worst','dP_orf_q50_worst','dP_orf_q95_worst', ...
              'T_oil_end_worst','T_steel_end_worst', ...
              'energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio','P_mech_sum' ...
            };
            % Provide aliases if needed
            try
                if ismember('T_end_worst', S.table.Properties.VariableNames) && ...
                   ~ismember('T_oil_end_worst', S.table.Properties.VariableNames)
                    S.table.T_oil_end_worst = S.table.T_end_worst; %#ok<AGROW>
                end
                if ismember('dP95_worst', S.table.Properties.VariableNames) && ...
                   ~ismember('dP_orf_q95_worst', S.table.Properties.VariableNames)
                    S.table.dP_orf_q95_worst = S.table.dP95_worst; %#ok<AGROW>
                end
                if ismember('E_orifice_sum', S.table.Properties.VariableNames) && ...
                   ismember('E_struct_sum', S.table.Properties.VariableNames) && ...
                   ~ismember('energy_tot_sum', S.table.Properties.VariableNames)
                    S.table.energy_tot_sum = S.table.E_orifice_sum + S.table.E_struct_sum; %#ok<AGROW>
                end
            catch
            end
            for jj = 1:numel(candCols)
                nm = candCols{jj};
                if ismember(nm, S.table.Properties.VariableNames)
                    val = S.table.(nm);
                    try
                        if ~all(isfinite(val(:)))
                            val(~isfinite(val)) = 0;
                        end
                    catch
                    end
                    meta.(nm) = val;
                end
            end
        end
    catch
    end

    memo(key) = meta;
    try, memo_store('set', key, meta); catch, end

    end

function xq = quant_clamp_x(x)
    % Apply the same quantization/clamps as in eval
    x = x(:)';
    x(1) = Utils.quantize_step(x(1),0.05);
    x(3) = Utils.quantize_step(x(3),0.05);
    x(4) = Utils.quantize_step(x(4),0.05);
    x(5) = Utils.quantize_step(x(5),0.05);
    x(6) = Utils.quantize_step(x(6),0.01);
    x(7) = Utils.quantize_step(x(7),0.02);
    x(2) = round(max(x(2),1));
    x(1) = min(max(x(1), 2.80), 3.60);
    x(2) = min(max(x(2), 5), 6);
    x(3) = min(max(x(3), 3.60), 4.00);
    x(4) = min(max(x(4), 3.80), 4.00);
    x(5) = min(max(x(5), 1.50), 3.60);
    x(6) = min(max(x(6), 0.95), 1.10);
    x(7) = min(max(x(7), 0.78), 0.90);
    xq = x;
end

function out = memo_store(cmd, key, val)
    persistent MEM;
    if isempty(MEM), MEM = containers.Map(); end
    switch lower(cmd)
        case 'get'
            if isKey(MEM,key), out = MEM(key); else, out = []; end
        case 'set'
            MEM(key) = val; out = true;
        case 'clear'
            MEM = containers.Map(); out = true;
        otherwise
            out = [];
    end
end

function P = decode_params_from_x(params0_, x_)
    d_o_mm = x_(1); n_orf = round(x_(2));
    g_lo = x_(3); g_mid = x_(4); g_hi = x_(5);
    PF_tau = x_(6); PF_gain = x_(7);
    P = params0_;
    P.orf.d_o = d_o_mm * 1e-3;         % mm -> m
    P.n_orf   = n_orf;
    P.A_o     = P.n_orf * (pi * P.orf.d_o^2 / 4);
    P.Qcap_big= max(P.orf.CdInf * P.A_o, 1e-9) * sqrt(2 * 1.0e9 / P.rho);
    n  = size(P.M,1); nStories = n-1;
    tg = ones(nStories,1) * g_mid;
    loN = min(3, nStories); if loN > 0, tg(1:loN) = g_lo; end
    hiN = min(2, nStories); if hiN > 0, tg(end-hiN+1:end) = g_hi; end
    P.toggle_gain = tg;
    if isfield(P,'cfg') && isfield(P.cfg,'PF')
        P.cfg.PF.tau  = PF_tau;
        P.cfg.PF.gain = PF_gain;
    end
end
