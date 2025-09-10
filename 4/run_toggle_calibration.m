function out = run_toggle_calibration(scaled, params, opts)
% RUN_TOGGLE_CALIBRATION: Tek değişkenli (toggle_gain) tarama ile en iyi g* seçimi.
%   OUT = RUN_TOGGLE_CALIBRATION(SCALED, PARAMS, OPTS)
%   - SCALED/PARAMS boşsa, parametreler ve kayıtlar otomatik yüklenir.
%   - OPTS.fields: g_min, g_max, g_step, mu_factors, mu_weights, thr

    if nargin < 3, opts = struct(); end

    % Veri ve parametreleri hazırla
    if nargin < 2 || isempty(params)
        try, setup; catch, end
        parametreler;
        params = struct('M',M,'C0',C0,'K',K,'k_sd',k_sd,'c_lam0',c_lam0,'Lori',Lori, ...
            'orf',orf,'rho',rho,'Ap',Ap,'A_o',A_o,'Qcap_big',Qcap_big,'mu_ref',mu_ref, ...
            'thermal',thermal,'T0_C',T0_C,'T_ref_C',T_ref_C,'b_mu',b_mu, ...
            'c_lam_min',c_lam_min,'c_lam_cap',c_lam_cap,'Lgap',Lgap, ...
            'cp_oil',cp_oil,'cp_steel',cp_steel,'steel_to_oil_mass_ratio',steel_to_oil_mass_ratio, ...
            'n_dampers_per_story',n_dampers_per_story,'toggle_gain',toggle_gain,'story_mask',story_mask, ...
            'resFactor',resFactor,'cfg',cfg,'n_orf',n_orf,'Dp',Dp,'d_w',d_w,'D_m',D_m,'n_turn',n_turn, ...
            'Kd',Kd,'Ebody',Ebody,'Gsh',Gsh,'story_height',story_height);
    end
    if nargin < 1 || isempty(scaled)
        [~, scaled] = load_ground_motions(Utils.try_warn(@() evalin('base','T1'),'T1 not in base')); %#ok<ASGLU>
    end

    % Ayarlar
    g_min = Utils.getfield_default(opts,'g_min',0.6);
    g_max = Utils.getfield_default(opts,'g_max',4);
    g_step= Utils.getfield_default(opts,'g_step',0.2);
    if isfield(opts,'g_list') && ~isempty(opts.g_list)
        g_list = opts.g_list(:).';
    elseif isfield(opts,'g') && ~isempty(opts.g)
        g_list = opts.g(:).';
    else
        g_list = g_min:g_step:g_max;
    end

    O = struct();
    O.mu_factors = Utils.getfield_default(opts,'mu_factors',[1.00]);
    O.mu_weights = Utils.getfield_default(opts,'mu_weights',[1]);
    O.thr = Utils.default_qc_thresholds(Utils.getfield_default(opts,'thr',struct()));
    O.thermal_reset = 'each'; O.order = 'natural'; O.quiet = true; O.do_export = false; O.use_orifice = true; O.use_thermal = true;

    % Tarama (paralel hızlandırma)
    ng = numel(g_list);
    f1m = nan(ng,1); f2m = nan(ng,1);
    dPviol = nan(ng,1); Qcapviol = nan(ng,1); cavviol = nan(ng,1); Tviol = nan(ng,1); muviol = nan(ng,1);
    Utils.try_warn(@() parpool_hard_reset(16), 'run_toggle_calibration: parpool');
    parfor i=1:ng
        Pi = params; Pi.toggle_gain = g_list(i);
        Si = run_batch_windowed(scaled, Pi, O);
        f1m(i) = mean(Si.table.PFA_w);
        f2m(i) = mean(Si.table.IDR_w);
        dPviol(i)   = mean(Si.table.dP95_worst   > O.thr.dP95_max);
        Qcapviol(i) = mean(Si.table.Qcap95_worst >= O.thr.Qcap95_max);
        cavviol(i)  = mean(Si.table.cav_pct_worst > max(O.thr.cav_pct_max,0));
        Tviol(i)    = mean(Si.table.T_end_worst  > O.thr.T_end_max);
        muviol(i)   = mean(Si.table.mu_end_worst < O.thr.mu_end_min);
    end
    v = struct('dP',dPviol,'Qcap',Qcapviol,'cav',cavviol,'T',Tviol,'mu',muviol);

    % Normalizasyon ve seçim (knee'ye yakınlık)
    f1n = (f1m - min(f1m)) ./ max(eps, (max(f1m)-min(f1m)));
    f2n = (f2m - min(f2m)) ./ max(eps, (max(f2m)-min(f2m)));
    d = hypot(f1n, f2n);
    [~,ix] = min(d);

    out = struct();
    out.g_list = g_list(:);
    out.f1_mean = f1m; out.f2_mean = f2m; out.d_norm = d;
    out.viol_rates = v;
    out.best_idx = ix; out.g_star = g_list(ix);
end
