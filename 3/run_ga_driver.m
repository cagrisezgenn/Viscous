function [X,F,gaout] = run_ga_driver(scaledOrSnap, params, optsEval, optsGA)
% === SAFE PARPOOL OPEN (cleanup + cap threads) ===
try
    parpool_hard_reset(16);
catch ME
    warning('[run_ga_driver] parpool init: %s', ME.message);
end
%RUN_GA_DRIVER Hybrid GA driver: accepts snapshot path or in-memory structs.
%   [X,F,GAOUT] = RUN_GA_DRIVER(SCALED_OR_PATH, PARAMS, OPTSEVAL, OPTSGA)
%   Two modes:
%     A) Snapshot: run_ga_driver('out/<ts>/snapshot.mat', params?, ...)
%        Loads 'scaled' (and params if not provided).
%     B) In-memory: run_ga_driver(scaled, params, ...)
%   Fitness evaluations have no IO/plots; no re-scaling is performed.

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
    if ~isfield(optsEval,'export'),     optsEval.export     = struct('plots',false,'ds',inf); end
    if ~isfield(optsEval,'debug'), optsEval.debug = true; end


    % ---------- GA options and objective ----------
    if nargin < 4 || isempty(optsGA), optsGA = struct; end
    rng(42);

    % decision vector: [d_o_mm, n_orf, g_lo, g_mid, g_hi, PF_tau, PF_gain]
    lb = [2.80, 5, 3.60, 3.80, 1.50, 0.95, 0.78];
    ub = [3.60, 6, 4.00, 4.00, 3.60, 1.10, 0.90];
    IntCon = 2;  % n_orf only

    % Pass GA population size into eval for debug windowing
    try, optsEval.popsize = Utils.getfield_default(optsGA,'PopulationSize',48); catch, end
    if ~isfield(optsEval,'debug') || isempty(optsEval.debug), optsEval.debug = false; end
    % Provide dataset signature for memo/debug keys
    try, dsig = sum([scaled.IM]) + sum([scaled.PGA]); catch, dsig = 0; end
    optsEval.dsig = dsig;
    obj = @(x) eval_design_fast(x, scaled, params, optsEval); % quantize/clamp inside

    % OutputFcn for generation-level debug reporting (works with parallel)
    outfun = @(options,state,flag) ga_output_debug(options,state,flag, params, optsEval);

    options = optimoptions('gamultiobj', ...
       'PopulationSize',    Utils.getfield_default(optsGA,'PopulationSize',56), ...
       'MaxGenerations',    Utils.getfield_default(optsGA,'MaxGenerations',50), ...
       'CrossoverFraction', Utils.getfield_default(optsGA,'CrossoverFraction',0.85), ...
       'MutationFcn',       Utils.getfield_default(optsGA,'MutationFcn',{@mutationadaptfeasible}), ...
       'ParetoFraction',    Utils.getfield_default(optsGA,'ParetoFraction',0.70), ...
       'StallGenLimit',     Utils.getfield_default(optsGA,'StallGenLimit',10), ...
       'DistanceMeasureFcn','distancecrowding', ...
       'UseParallel',       Utils.getfield_default(optsGA,'UseParallel',true), ...
       'Display','iter','PlotFcn',[], 'FunctionTolerance',1e-3, ...
       'OutputFcn', outfun);

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
    outdir = fullfile('out', ['ga_' datestr(now,'yyyy-mm-dd_HHMMSS')]);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    opts_ga = options; date_str = datestr(now);
    save(fullfile(outdir,'ga_front.mat'),'X','F','opts_ga','meta','date_str','-v7.3');

    % CSV (design vars + objectives + penalty parts from memo)
    pen     = nan(size(X,1),1);
    pen_dP  = nan(size(X,1),1);
    pen_Qcap= nan(size(X,1),1);
    pen_cav = nan(size(X,1),1);
    pen_T   = nan(size(X,1),1);
    pen_mu  = nan(size(X,1),1);
    try
        for i = 1:size(X,1)
            Xi = quant_clamp_x(X(i,:));
            try, dsig_pack = sum([scaled.IM]) + sum([scaled.PGA]); catch, dsig_pack = 0; end
            key_i = jsonencode([Xi, dsig_pack]);
            mi = memo_store('get', key_i);
            if ~isempty(mi)
                pen(i,1) = Utils.getfield_default(mi,'pen',NaN);
                pp = Utils.getfield_default(mi,'pen_parts',struct());
                pen_dP(i,1)   = Utils.getfield_default(pp,'dP',NaN);
                pen_Qcap(i,1) = Utils.getfield_default(pp,'Qcap',NaN);
                pen_cav(i,1)  = Utils.getfield_default(pp,'cav',NaN);
                pen_T(i,1)    = Utils.getfield_default(pp,'T',NaN);
                pen_mu(i,1)   = Utils.getfield_default(pp,'mu',NaN);
            end
        end
    catch
    end
    T = array2table([X F pen pen_dP pen_Qcap pen_cav pen_T pen_mu], 'VariableNames', ...
        {'d_o_mm','n_orf','g_lo','g_mid','g_hi','PF_tau','PF_gain','f1','f2', ...
         'pen','pen_dP','pen_Qcap','pen_cav','pen_T','pen_mu'});
    writetable(T, fullfile(outdir,'ga_front.csv'));

    % Knee-point (min distance to ideal after min-max normalization)
    try
        Tfront = readtable(fullfile(outdir,'ga_front.csv'));
        f1v = Tfront.f1; f2v = Tfront.f2;
        f1n = (f1v - min(f1v)) / max(eps, (max(f1v)-min(f1v)));
        f2n = (f2v - min(f2v)) / max(eps, (max(f2v)-min(f2v)));
        d   = hypot(f1n, f2n);
        [~,kidx] = min(d);
        Tknee = Tfront(kidx,:);
        writetable(Tknee, fullfile(outdir,'ga_knee.csv'));
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

    % Minimal README
    fid=fopen(fullfile(outdir,'README.txt'),'w');
    if fid~=-1
        fprintf(fid, 'Hybrid GA run: %s\n', date_str);
        fprintf(fid, 'IM_mode=%s, band_fac=%s, s_bounds=%s\n', ...
            mat2str(meta.IM_mode), mat2str(meta.band_fac), mat2str(meta.s_bounds));
        fprintf(fid, 'mu_factors=%s, mu_weights=%s\n', mat2str(meta.mu_factors), mat2str(meta.mu_weights));
        try
            fprintf(fid, 'thr=%s\n', jsonencode(meta.thr));
        catch
        end
        fprintf(fid, 'Note: No simulations during packaging. Fitness evals had no IO.\n');
        fclose(fid);
    end
end

function [state, options, optchanged] = ga_output_debug(options, state, flag, params, optsEval)
%GA_OUTPUT_DEBUG OutputFcn to report unique designs per generation.
    optchanged = false;
    try
        if ~isfield(state,'Population') || isempty(state.Population)
            return;
        end
        if ~(strcmpi(flag,'iter') || strcmpi(flag,'done'))
            return;
        end
        X = state.Population;
        n = size(X,1);
        Xq = nan(n,7);
        for i = 1:n
            Xq(i,:) = quant_clamp_x(X(i,:));
        end
        % derive params signatures: [xq, A_o, Qcap_big]
        d_o_m = Xq(:,1) * 1e-3;
        n_orf = round(Xq(:,2));
        Ao = n_orf .* (pi * (d_o_m.^2) / 4);
        Qcap_big = max(params.orf.CdInf .* Ao, 1e-9) .* sqrt(2 * 1.0e9 / params.rho);
        Pmat = [Xq Ao Qcap_big];
        % unique counts
        uX = unique(Xq, 'rows');
        uP = unique(Pmat, 'rows');
        ucnt_x = size(uX,1);
        ucnt_p = size(uP,1);
        popsz = n;
        % Diversity threshold and penalty share
        thr = 10; if isstruct(optsEval) && isfield(optsEval,'uniq_thr') && ~isempty(optsEval.uniq_thr), thr = optsEval.uniq_thr; end
        dsig = 0; if isstruct(optsEval) && isfield(optsEval,'dsig'), dsig = optsEval.dsig; end
        share_thr = 0.5; if isstruct(optsEval) && isfield(optsEval,'pen_share_thr') && ~isempty(optsEval.pen_share_thr), share_thr = optsEval.pen_share_thr; end
        pen_hi = 0;
        for i = 1:n
            key_i = jsonencode([Xq(i,:), dsig]);
            mi = memo_store('get', key_i);
            if ~isempty(mi)
                share_i = Utils.getfield_default(mi,'pen',NaN);
                if isfinite(share_i) && (share_i > share_thr)
                    pen_hi = pen_hi + 1;
                end
            end
        end
        pen_rate = pen_hi / max(popsz,1);
        flagLow = (ucnt_p < thr);
        fprintf('[GA dbg/of] gen=%d | unique x=%d/%d, unique params=%d/%d | high-pen-rate=%.0f%% (thr=%.2f)%s\n', ...
            state.Generation, ucnt_x, popsz, ucnt_p, popsz, 100*pen_rate, share_thr, tern(flagLow,'  <-- LOW DIVERSITY',''));
    catch
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
    O.export = struct('plots', false, 'ds', inf);
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
        return;
    end

    % --- HARD FILTER (erken eleme) ---
    dP95v = S.table.dP95_worst;
    qcapv = S.table.Qcap95_worst;
    cavv  = S.table.cav_pct_worst;
    if any(dP95v > 1e9) || any(qcapv > 0.90) || any(cavv > 0.01)
        f = [1e6, 1e6];
        meta = struct('x',x,'f',f,'hard_kill',true);
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
    meta = struct('x',x,'f',f,'PFA_w_mean',f1,'IDR_w_mean',f2,'pen',lambda*pen, ...
                  'pen_parts',struct('dP',pen_dP,'Qcap',pen_Qcap,'cav',pen_cav,'T',pen_T,'mu',pen_mu));

    % ---------------- Debug/pilot logging: uniqueness and penalty share ----------------
    try
        dbg_enabled = true;
        if nargin >= 4 && ~isempty(optsEval) && isstruct(optsEval)
            dbg_enabled = isfield(optsEval,'debug') && logical(optsEval.debug);
        end
        if dbg_enabled
            % Build signatures:
            %  - x_sig: quantized/clamped decision vector (x)
            %  - p_sig: used param tuple [x, P.A_o, P.Qcap_big]
            x_sig = mat2str(x, 6);
            try
                p_sig_vec = [x(1:7) P.A_o P.Qcap_big];
            catch
                p_sig_vec = [x(1:7) NaN NaN];
            end
            p_sig = mat2str(p_sig_vec, 6);
            % Persistent state across calls
            persistent DBG;
            if isempty(DBG)
                DBG = struct();
                DBG.popsize = Utils.getfield_default(optsEval,'popsize',24);
                DBG.count = 0;                % total eval calls since last reset
                DBG.sigset_x = containers.Map(); % signatures (x) in current generation window
                DBG.sigset_p = containers.Map(); % signatures (params) in current generation window
                DBG.pen_hi = 0;               % high-penalty counter in window
                DBG.share_thr = Utils.getfield_default(optsEval,'pen_share_thr',0.5);
            end
            % Start of a generation window?
            if mod(DBG.count, DBG.popsize) == 0
                DBG.sigset_x = containers.Map();
                DBG.sigset_p = containers.Map();
                DBG.pen_hi = 0;
            end
            % Add signature
            if ~isKey(DBG.sigset_x, x_sig)
                DBG.sigset_x(x_sig) = 1;
            else
                DBG.sigset_x(x_sig) = DBG.sigset_x(x_sig) + 1;
            end
            if ~isKey(DBG.sigset_p, p_sig)
                DBG.sigset_p(p_sig) = 1;
            else
                DBG.sigset_p(p_sig) = DBG.sigset_p(p_sig) + 1;
            end
            % Penalty share (use multiplicative share lambda*pen)
            share = lambda * pen;
            if (share > DBG.share_thr)
                DBG.pen_hi = DBG.pen_hi + 1;
            end
            DBG.count = DBG.count + 1;
            % End of generation window -> print summary
            if mod(DBG.count, DBG.popsize) == 0
                genIdx = DBG.count / DBG.popsize;
                try, uniq_x = DBG.sigset_x.Count; catch, uniq_x = numel(keys(DBG.sigset_x)); end
                try, uniq_p = DBG.sigset_p.Count; catch, uniq_p = numel(keys(DBG.sigset_p)); end
                pen_rate = DBG.pen_hi / max(DBG.popsize,1);
                fprintf('[GA dbg] gen=%d | unique x=%d/%d, unique params=%d/%d | high-pen-rate=%.0f%% (thr=%.2f)\n', ...
                    genIdx, uniq_x, DBG.popsize, uniq_p, DBG.popsize, 100*pen_rate, DBG.share_thr);
            end
        end
    catch
    end
    memo(key) = meta;

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
