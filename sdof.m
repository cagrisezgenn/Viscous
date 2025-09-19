%% SDOF analysis using GA best parameters
% Runs a single simulation without optimization using the best design
% parameters (xbest) stored in ga_front.csv. Produces displacement and
% acceleration histories for the 10th floor along with peak metrics.

clear; clc; close all;

%% Load base parameters (embedded defaults)
params = default_params();

M = params.M;        C0 = params.C0;       K = params.K;
n = params.n;        story_height = params.story_height;
T1 = params.T1;

Dp = params.Dp;      Lgap = params.Lgap;   d_w = params.d_w;
D_m = params.D_m;    n_turn = params.n_turn;
mu_ref = params.mu_ref;  Lori = params.Lori;
n_orf = params.n_orf;    rho = params.rho;
orf = params.orf;        cfg = params.cfg;
thermal = params.thermal;
T0_C = params.T0_C;  T_ref_C = params.T_ref_C;  b_mu = params.b_mu;
cp_oil = params.cp_oil;   cp_steel = params.cp_steel;
steel_to_oil_mass_ratio = params.steel_to_oil_mass_ratio;
n_dampers_per_story = params.n_dampers_per_story;
story_mask = params.story_mask;  resFactor = params.resFactor;
Kd = params.Kd; Ebody = params.Ebody; Gsh = params.Gsh;

Ap = params.Ap;   Ao = params.Ao;   k_sd = params.k_sd;
c_lam0 = params.c_lam0;  c_lam_cap = params.c_lam_cap;
c_lam_min_abs = params.c_lam_min_abs;  c_lam_min_frac = params.c_lam_min_frac;
c_lam_min = params.c_lam_min;
Qcap_big = params.Qcap_big;

%% Override with GA best parameters from GA output
try
    tbl = readtable('ga_front.csv');
    xb = tbl(1,:);

    d_o = xb.d_o_mm/1000;      % [m]
    n_orf = xb.n_orf;
    orf.Cd0  = xb.Cd0;
    orf.CdInf= xb.CdInf;
    orf.p_exp= xb.p_exp;
    Lori  = xb.Lori_mm/1000;    % [m]
    thermal.hA_W_perK = xb.hA_W_perK;
    Dp    = xb.Dp_mm/1000;      % [m]
    d_w   = xb.d_w_mm/1000;     % [m]
    D_m   = xb.D_m_mm/1000;     % [m]
    n_turn= xb.n_turn;
    mu_ref= xb.mu_ref;

    orf.d_o = d_o;

    params_override = struct('Dp',Dp,'d_w',d_w,'D_m',D_m,'n_turn',n_turn, ...
        'mu_ref',mu_ref,'Lori',Lori,'Lgap',Lgap,'Kd',Kd,'Ebody',Ebody, ...
        'Gsh',Gsh,'orf',orf,'n_orf',n_orf,'rho',rho, ...
        'c_lam_min_abs',c_lam_min_abs,'c_lam_min_frac',c_lam_min_frac);

    params_override.cfg = cfg;
    params_override = build_params_local(params_override);

    Ap = params_override.Ap;
    Ao = params_override.Ao;
    k_sd = params_override.k_sd;
    c_lam0 = params_override.c_lam0;
    c_lam_min = params_override.c_lam_min;
    Qcap_big = params_override.Qcap_big;
    orf = params_override.orf;
catch ME
    warning('Unable to apply GA parameters: %s', ME.message);
end

%% Ground motion input (scaled by T1)
[~, recs] = load_ground_motions_local(T1);
rec = recs(1);  t = rec.t; ag = rec.ag;
win = make_arias_window_local(t, ag); t5 = win.t5; t95 = win.t95;

%% Damperless and damper responses
[x0,a_rel0] = solve_linear_mck_local(t, ag, M, C0, K);
[x_d,a_d,diag_d] = mck_with_damper_local(t, ag, M, C0, K, k_sd, c_lam0, Lori, ...
    orf, rho, Ap, Ao, Qcap_big, mu_ref, thermal, T0_C, T_ref_C, b_mu, ...
    c_lam_min, c_lam_cap, Lgap, cp_oil, cp_steel, steel_to_oil_mass_ratio, ...
    story_mask, n_dampers_per_story, resFactor, cfg);

x10_0 = x0(:,10);
x10_d = x_d(:,10);
a10_0 = a_rel0(:,10) + ag;
a10_d = a_d(:,10) + ag;

%% Assemble equivalent damping/stiffness matrices for modal check
% nStories = n - 1;
% mask = story_mask(:); if numel(mask)==1, mask = mask*ones(nStories,1); end
% ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
% multi = (mask .* ndps);
% Kadd = zeros(n); C_add = zeros(n);
% for i = 1:nStories
%     idx = [i,i+1];
%     k_eq = k_sd * multi(i);
%     c_eq = diag_d.c_lam * multi(i);
%     kM  = k_eq * [1 -1; -1 1];
%     cM  = c_eq * [1 -1; -1 1];
%     Kadd(idx,idx) = Kadd(idx,idx) + kM;
%     C_add(idx,idx)= C_add(idx,idx) + cM;
% end
% K_tot = K + Kadd;
% C_d = C0 + C_add;
% [V,D] = eig(K_tot,M); [w2,ord] = sort(diag(D),'ascend');
% phi1 = V(:,ord(1)); w1 = sqrt(w2(1));
% normM = phi1.' * M * phi1;

[V_base,D_base] = eig(K,M); [w2_base,ord_base] = sort(diag(D_base),'ascend');
phi1 = V_base(:,ord_base(1));
w1 = sqrt(w2_base(1));
normM = phi1.' * M * phi1;

%% Plots: 10th floor displacement and acceleration, plus IDR
figure('Name','10. Kat yer de??i??tirme ??? ham ivme (ODE-only)','Color','w');
plot(t, x10_0,'k','LineWidth',1.4); hold on;
plot(t, x10_d,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('x10(t) [m]');
title(sprintf('10-Kat | T1=%.3f s | Arias [%.3f, %.3f] s', T1, t5, t95));
legend('Dampersiz','Damperli','Location','best');

figure('Name','10. Kat mutlak ivme','Color','w');
plot(t, a10_0,'k','LineWidth',1.4); hold on;
plot(t, a10_d,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('a10abs(t) [m/s^2]');
legend('Dampersiz','Damperli','Location','best');

drift0 = x0(:,2:end) - x0(:,1:end-1);
drift_d = x_d(:,2:end) - x_d(:,1:end-1);
idx = win.idx;
IDR0 = max(abs(drift0(idx,:)))./story_height;
IDR_d = max(abs(drift_d(idx,:)))./story_height;
story_ids = 1:(n-1);
figure('Name','Maksimum IDR','Color','w');
plot(story_ids, IDR0,'k-o','LineWidth',1.4); hold on;
plot(story_ids, IDR_d,'r-d','LineWidth',1.0);
grid on; xlabel('Kat'); ylabel('Maks IDR [Delta x/h]');
legend('Dampersiz','Damperli','Location','best');

%% Summary printout
zeta0 = (phi1.' * C0 * phi1) / (2*w1*normM);

idx = win.idx;
x10_max_0 = max(abs(x10_0(idx)));
x10_max_damperli = max(abs(x10_d(idx)));
a10abs_max_0 = max(abs(a10_0(idx)));
a10abs_max_damperli = max(abs(a10_d(idx)));
IDR_max_0 = max(IDR0);
IDR_max_damperli = max(IDR_d);

fprintf(['Self-check zeta1 (yapısal baz): %.3f%%; ' ...
         'x10_{max}0=%.4g m; x10_{max}d=%.4g m; ' ...
         'a10abs_{max}0=%.4g m/s^2; a10abs_{max}d=%.4g m/s^2; ' ...
         'IDR_{max}0=%.4g; IDR_{max}d=%.4g\n'], ...
        100*zeta0, x10_max_0, x10_max_damperli, ...
        a10abs_max_0, a10abs_max_damperli, IDR_max_0, IDR_max_damperli);

%% Metric aggregation and CSV export  (CLEAN)
metrics = compute_sdof_metrics_local(t, x_d, a_d, ag, diag_d, story_height, win, Qcap_big);

% SCI-rapor metrik adları (mean alanlarını aynı değere eşle)
metrics.PFA_mean = metrics.PFA;
metrics.IDR_mean = metrics.IDR;

% Ceza terimleri
metrics = compute_penalties_local(metrics);

% Artık dP95/dP50, Q_q50, Q_q95, Qcap95, PF_p95, cav_pct, T_end, mu_end
% hepsi compute_sdof_metrics_local içinden geliyor — tekrar hesaplama YOK.

row = struct( ...
    'PFA_mean', metrics.PFA_mean, ...
    'IDR_mean', metrics.IDR_mean, ...
    'pen', metrics.pen, ...
    'pen_dP', metrics.pen_dP, ...
    'pen_Qcap', metrics.pen_Qcap, ...
    'pen_cav', metrics.pen_cav, ...
    'pen_T', metrics.pen_T, ...
    'pen_mu', metrics.pen_mu, ...
    'x10_max_damperli', metrics.x10_max_damperli, ...
    'a10abs_max_damperli', metrics.a10abs_max_damperli, ...
    'dP95', metrics.dP95, ...
    'Qcap95', metrics.Qcap95, ...
    'cav_pct', metrics.cav_pct, ...
    'T_end', metrics.T_end, ...
    'mu_end', metrics.mu_end, ...
    'PF_p95', metrics.PF_p95, ...
    'Q_q50', metrics.Q_q50, ...
    'Q_q95', metrics.Q_q95, ...
    'dP50', metrics.dP50, ...
    'energy_tot_sum', metrics.energy_tot_sum, ...
    'E_orifice_sum', metrics.E_orifice_sum, ...
    'E_struct_sum', metrics.E_struct_sum, ...
    'E_ratio', metrics.E_ratio, ...
    'P_mech_sum', metrics.P_mech_sum );

T_metrics = struct2table(row);
if ~exist('out','dir'); mkdir('out'); end
tstamp  = datestr(now,'yyyymmdd_HHMMSS');
outfile = fullfile('out', ['sdof_' tstamp '.csv']);
writetable(T_metrics, outfile);
fprintf('SDOF metrics written to %s\n', outfile);
%% Local helpers
function metrics = compute_sdof_metrics_local(t, x_d, a_rel_d, ag, diag_d, story_height, win, Qcap_big)
    if nargin < 8 || ~isfinite(Qcap_big) || Qcap_big <= 0
        Qcap_big = eps;
    end
    idx = win.idx(:);
    if isempty(idx) || all(~idx)
        idx = true(size(t));
    end
    a_top_abs = a_rel_d(:,end) + ag(:);
    metrics = struct();
    metrics.PFA = max(abs(a_top_abs(idx)));
    metrics.x10_max_damperli = max(abs(x_d(idx,end)));
    metrics.a10abs_max_damperli = metrics.PFA;
    drift = (x_d(:,2:end) - x_d(:,1:end-1)) / story_height;
    metrics.IDR = max(max(abs(drift(idx,:))));
    if ~isfinite(metrics.IDR)
        metrics.IDR = 0;
    end
    abs_dP = abs(diag_d.dP_eff(idx,:));
    abs_story_force = abs(diag_d.story_force(idx,:));
    dP_q50 = quantile(abs_dP, 0.50);
    dP_q95 = quantile(abs_dP, 0.95);
    story_force_q95 = quantile(abs_story_force, 0.95);
    [~, ws] = max(story_force_q95);
    if isempty(ws) || ~isfinite(ws)
        ws = 1;
    end
    metrics.dP95 = dP_q95(ws);
    metrics.dP50 = dP_q50(ws);
    if isfield(diag_d,'Q') && ~isempty(diag_d.Q)
        abs_Q = abs(diag_d.Q(idx,:));
        Q_q50 = quantile(abs_Q, 0.50);
        Q_q95 = quantile(abs_Q, 0.95);
        metrics.Q_q50 = Q_q50(ws);
        metrics.Q_q95 = Q_q95(ws);

        Qcap_ratio = abs_Q ./ max(Qcap_big, eps);
        Qcap95_vec = quantile(Qcap_ratio, 0.95);
        metrics.Qcap95 = Qcap95_vec(ws);
    else
        [metrics.Q_q50, metrics.Q_q95, metrics.Qcap95] = deal(0);
    end
    cav_mean = mean(diag_d.cav_mask(idx,:), 1);
    metrics.cav_pct = cav_mean(ws);
    E_orf_end    = isfield(diag_d,'E_orf')    && ~isempty(diag_d.E_orf);
    E_struct_end = isfield(diag_d,'E_struct') && ~isempty(diag_d.E_struct);
    if E_orf_end
        metrics.E_orifice_sum = diag_d.E_orf(end);
    else
        metrics.E_orifice_sum = 0;
    end
    if E_struct_end
        metrics.E_struct_sum = diag_d.E_struct(end);
    else
        metrics.E_struct_sum = 0;
    end
    metrics.energy_tot_sum = metrics.E_orifice_sum + metrics.E_struct_sum;
    if metrics.E_struct_sum > 0
        metrics.E_ratio = metrics.E_orifice_sum / max(metrics.E_struct_sum, eps);
    else
        metrics.E_ratio = 0;
    end
    if isfield(diag_d,'P_sum') && ~isempty(diag_d.P_sum)
        P_mech_win = diag_d.P_sum(idx);
        P_mech_win = P_mech_win(isfinite(P_mech_win));
        metrics.P_mech_sum = mean(P_mech_win);
    else
        metrics.P_mech_sum = NaN;
    end
    if isfield(diag_d,'PF') && ~isempty(diag_d.PF)
        PF_abs = abs(diag_d.PF(idx,:));
        if size(PF_abs,2) >= ws
            metrics.PF_p95 = quantile(PF_abs(:,ws), 0.95);
        else
            metrics.PF_p95 = NaN;
        end
    else
        metrics.PF_p95 = NaN;
    end
    if isfield(diag_d,'T_oil') && ~isempty(diag_d.T_oil)
        metrics.T_end = diag_d.T_oil(end);
    else
        metrics.T_end = NaN;
    end
    if isfield(diag_d,'mu') && ~isempty(diag_d.mu)
        metrics.mu_end = diag_d.mu(end);
    else
        metrics.mu_end = NaN;
    end

end
function y = softmin_local(a,b,epsm)
    if nargin < 3 || isempty(epsm)
        epsm = 1e5;
    end
    y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
end

function win = make_arias_window_local(t, ag, varargin)
    p = inputParser;
    p.addParameter('p1',0.05,@(x)isscalar(x) && x>=0 && x<=1);
    p.addParameter('p2',0.95,@(x)isscalar(x) && x>=0 && x<=1);
    p.addParameter('pad',0.5,@(x)isscalar(x) && x>=0);
    p.parse(varargin{:});
    p1 = p.Results.p1; p2 = p.Results.p2; pad = p.Results.pad;

    IA = cumtrapz(t, ag.^2);
    IA_tot = IA(end);
    if ~(isfinite(IA_tot)) || IA_tot <= eps
        t_start = t(1); t_end = t(end);
        idx = true(size(t));
        win = struct('t5',t_start,'t95',t_end,'pad',0, ...
                     't_start',t_start,'t_end',t_end,'idx',idx, ...
                     'coverage',1.0,'flag_low_arias',true);
        return;
    end

    IA_norm = IA / IA_tot;
    [t_unique, iu] = unique(IA_norm);
    interp_t = t(iu);
    t5  = interp1(t_unique, interp_t, p1, 'linear');
    t95 = interp1(t_unique, interp_t, p2, 'linear');
    dur = t95 - t5;
    if dur < 5, pad = 0.25; end
    t_start = max(t(1),  t5  - pad);
    t_end   = min(t(end), t95 + pad);
    idx = (t >= t_start) & (t <= t_end);
    coverage = trapz(t(idx), ag(idx).^2) / IA_tot;
    flag_low_arias = coverage < 0.90;
    win = struct('t5',t5,'t95',t95,'pad',pad, ...
                 't_start',t_start,'t_end',t_end,'idx',idx, ...
                 'coverage',coverage,'flag_low_arias',flag_low_arias);
end

function [records, scaled, meta] = load_ground_motions_local(T1, opts)
    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'verbose'), opts.verbose = true; end
    hp_cut   = getfield_default_local(opts,'hp_cut',0.05);
    IM_mode  = getfield_default_local(opts,'IM_mode','band');
    band_fac = getfield_default_local(opts,'band_fac',[0.8 1.2]);
    band_N   = getfield_default_local(opts,'band_N',21);
    s_bounds = getfield_default_local(opts,'s_bounds',[0.2 2.2]);

    if exist('acc_matrix.mat','file')
        raw = load('acc_matrix.mat');
    else
        warning('load_ground_motions_local:synthetic', ...
            'acc_matrix.mat bulunamad??; sentetik tek kay??t kullan??lacak.');
        t_syn = (0:0.01:40)';
        ag_syn = 0.35 * sin(2*pi*1.2*t_syn) .* exp(-0.05*t_syn);
        raw = struct('synthetic',[t_syn ag_syn]);
    end
    fn  = fieldnames(raw);
    records = struct('name',{},'t',{},'ag',{},'dt',{},'duration',{},'PGA',{},'PGV',{},'IM',{},'scale',{});

    for kk = 1:numel(fn)
        A = raw.(fn{kk});
        t  = A(:,1);   ag = A(:,2);
        [t,iu] = unique(t,'stable'); ag = ag(iu);
        dt = median(diff(t));
        assert(max(abs(diff(t) - dt)) < 1e-6, 'Zaman ??rnekleme aral?????? d??zensiz');
        t  = (t(1):dt:t(end)).';
        ag = interp1(A(:,1), A(:,2), t, 'linear');
        assert(max(abs(ag)) < 100, '??vme b??y??kl?????? birim hatas??na i??aret ediyor');
        ag = detrend(ag,0);
        ag = detrend(ag,1);
        try
            Fs = 1/dt;
            Wn = hp_cut/(Fs/2);
            [b,a] = butter(2,Wn,'high');
            ag = filtfilt(b,a,ag);
        catch
        end
        duration = t(end) - t(1);
        v  = cumtrapz(t,ag);
        PGA = max(abs(ag));
        PGV = max(abs(v));
        records(end+1) = struct('name',fn{kk},'t',t,'ag',ag,'dt',dt, ...
                                 'duration',duration,'PGA',PGA,'PGV',PGV, ...
                                 'IM',[],'scale',1); %#ok<AGROW>
    end

    if getfield_default_local(opts,'verbose',true)
        fprintf('Toplam %d zemin hareketi kayd?? y??klendi:\n', numel(records));
        for kk = 1:numel(records)
            r = records(kk);
            fprintf('%2d) %-12s dt=%6.4f s dur=%6.2f s PGA=%7.3f PGV=%7.3f\n', ...
                kk, r.name, r.dt, r.duration, r.PGA, r.PGV);
        end
    end

    scaled = [];
    meta = struct();

    if nargin >= 1 && ~isempty(T1)
        for kk = 1:numel(records)
            records(kk).IM = compute_IM_local(records(kk).t, records(kk).ag, IM_mode, T1, band_fac, band_N);
        end

        IM = [records.IM];
        IM_low  = max(s_bounds(1)*IM);
        IM_high = min(s_bounds(2)*IM);
        targetIM0 = median(IM);
        targetIM  = min(max(targetIM0, IM_low), IM_high);
        doClip = (IM_low > IM_high);

        max_iter = 3; dropped = {};
        for it = 1:max_iter
            IM_low  = max(s_bounds(1)*IM);
            IM_high = min(s_bounds(2)*IM);
            if IM_low <= IM_high, break; end
            [~,idx] = max(abs(log(IM) - median(log(IM))));
            dropped{end+1} = records(idx).name; %#ok<AGROW>
            records(idx) = [];
            IM(idx) = [];
        end
        if ~isempty(dropped) && getfield_default_local(opts,'verbose',true)
            fprintf('TRIM: ay??klanan u?? de??erler = %s\n', strjoin(dropped,', '));
        end
        IM_low  = max(s_bounds(1)*IM);
        IM_high = min(s_bounds(2)*IM);
        targetIM0 = median(IM);
        targetIM  = min(max(targetIM0, IM_low), IM_high);
        doClip = (IM_low > IM_high);

        scaled = records; n_clipped = 0; s_all = zeros(1,numel(records));
        for kk = 1:numel(records)
            s_raw = targetIM / records(kk).IM;
            if doClip
                s = min(max(s_raw, s_bounds(1)), s_bounds(2));
                if abs(s - s_raw) > 1e-12, n_clipped = n_clipped + 1; end
            else
                s = s_raw;
            end
            s_all(kk) = s;

            scaled(kk).ag    = s * records(kk).ag;
            scaled(kk).PGA   = s * records(kk).PGA;
            scaled(kk).PGV   = s * records(kk).PGV;
            scaled(kk).scale = s;
            scaled(kk).s_clipped = doClip && (abs(s - s_raw) > 1e-12);
            scaled(kk).trimmed   = false;
            scaled(kk).IM = compute_IM_local(scaled(kk).t, scaled(kk).ag, IM_mode, T1, band_fac, band_N);
        end

        err = abs([scaled.IM] - targetIM) / max(targetIM, eps) * 100;
        modeStr = tern_local(strcmpi(IM_mode,'band'),'band','PSA@T1');
        clipCount = n_clipped * doClip;
        if getfield_default_local(opts,'verbose',true)
            fprintf('Hedef IM = %.3f (%s). Maks hata = %.2f%% | uygun aral??k=[%.3f, %.3f] | s_min=%.2f s_max=%.2f | KIRPILAN=%d\n', ...
                targetIM, modeStr, max(err), IM_low, IM_high, min(s_all), max(s_all), clipCount);
        end
        meta = struct('IM_mode', IM_mode, 'band_fac', band_fac, 's_bounds', s_bounds);
        if exist('dropped','var'), meta.TRIM_names = dropped; else, meta.TRIM_names = {}; end
    end
end

function IM = compute_IM_local(t, ag, mode, T1, band_fac, band_N)
    zeta = 0.05;
    if strcmpi(mode,'band')
        Tgrid = linspace(band_fac(1)*T1, band_fac(2)*T1, band_N);
        Sa = zeros(size(Tgrid));
        for ii = 1:numel(Tgrid)
            Sa(ii) = calc_psa_local(t, ag, Tgrid(ii), zeta);
        end
        IM = exp(mean(log(Sa + eps)));
    else
        IM = calc_psa_local(t, ag, T1, zeta);
    end
end

function Sa = calc_psa_local(t, ag, T, zeta)
    w = 2*pi / T;
    agf = griddedInterpolant(t, ag, 'linear', 'nearest');
    odef = @(tt, y)[ y(2);
                     -2*zeta*w*y(2) - w*w*y(1) - agf(tt) ];
    y0 = [0;0];
    [~, y] = ode45(odef, t, y0);
    x  = y(:,1);
    xd = y(:,2);
    xdd = -2*zeta*w*xd - w*w*x - ag;
    abs_acc = xdd + ag;
    Sa = max(abs(abs_acc));
end

function s = tern_local(c,a,b)
    if c, s=a; else, s=b; end
end

function [x,a_rel] = solve_linear_mck_local(t, ag, M, C, K)
    nloc = size(M,1);
    r = ones(nloc,1);
    ag_interp = griddedInterpolant(t, ag, 'linear', 'nearest');
    odefun = @(tt,z) [ z(nloc+1:end); ...
        M \ (-C*z(nloc+1:end) - K*z(1:nloc) - M*r*ag_interp(tt)) ...
    ];
    z0 = zeros(2*nloc,1);
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    sol = ode15s(odefun, [t(1) t(end)], z0, opts);
    z = deval(sol, t).';
    x = z(:,1:nloc);
    a_rel = (-(M \ (C*z(:,nloc+1:end).' + K*z(:,1:nloc).')).' - ag.*r.');
end

function [x,a_rel,ts] = mck_with_damper_local(t,ag,M,C,K, k_sd,c_lam0,Lori, orf,rho,Ap,Ao,Qcap, mu_ref, ...
    thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, story_mask, ...
    n_dampers_per_story, resFactor, cfg)

    Nt = numel(t);
    n = size(M,1);
    if n < 2
        x = zeros(Nt,n);
        a_rel = zeros(Nt,n);
        ts = struct();
        return;
    end

    Ns = n - 1;
    r = ones(n,1);
    agf = griddedInterpolant(t, ag, 'linear', 'nearest');

    mask = story_mask(:);
    if numel(mask)==1, mask = mask * ones(Ns,1); end
    ndps = n_dampers_per_story(:);
    if numel(ndps)==1, ndps = ndps * ones(Ns,1); end
    multi_col = max(mask .* ndps, 0);
    multi_row = multi_col.';
    Nvec = 1:Ns;
    Mvec = 2:n;

    cfg_on = getfield_default_local(cfg,'on', struct());
    mu_floor_flag = getfield_default_local(cfg_on,'mu_floor', false);
    pf_res_only = getfield_default_local(cfg_on,'pf_resistive_only', false);

    numStruct = getfield_default_local(cfg,'num', struct());
    mu_min_phys = getfield_default_local(numStruct,'mu_min_phys', mu_ref);
    mu_max_phys = getfield_default_local(numStruct,'mu_max_phys', NaN);
    soft_eps = getfield_default_local(numStruct,'softmin_eps', 1e5);
    Qcap_scale = getfield_default_local(numStruct,'Qcap_scale', 1);

    c_lam = min(max(c_lam0, c_lam_min), c_lam_cap);

    PF_cfg = getfield_default_local(cfg,'PF', struct());
    gain_pf = getfield_default_local(PF_cfg,'gain',1);

    rho_base = getfield_default_local(thermal,'rho_ref', rho);
    alpha_rho = getfield_default_local(thermal,'alpha_rho', 0);
    beta0 = getfield_default_local(thermal,'beta0', 1.6e9);
    b_beta = getfield_default_local(thermal,'b_beta', 0);
    T_ref = getfield_default_local(thermal,'T_ref_C', T_ref_C);
    T_env = getfield_default_local(thermal,'T_env_C', T0_C);
    dT_max = getfield_default_local(thermal,'dT_max', 80);

    Ap_single = Ap;
    Ao_single = Ao;
    Ap_story_col = Ap_single * multi_col;
    Ao_story_col = Ao_single * multi_col;
    Ap_story_row = Ap_story_col.';
    Ao_story_row = Ao_story_col.';
    n_orf_single = getfield_default_local(orf,'n_orf',1);
    n_orf_total_col = max(1, n_orf_single) * multi_col;
    n_orf_total_row = n_orf_total_col.';
    d_o_single = getfield_default_local(orf,'d_o',3.0e-3);

    Qcap_eff = max(1e-9, Qcap * Qcap_scale);
    Qcap_story = Qcap_eff * ones(1,Ns);

    resFactor_loc = getfield_default_local(thermal,'resFactor', resFactor);
    cp_oil_loc = getfield_default_local(thermal,'cp_oil', cp_oil);
    cp_steel_loc = getfield_default_local(thermal,'cp_steel', cp_steel);
    hA_base = getfield_default_local(thermal,'hA_W_perK', thermal.hA_W_perK);
    hA_os = getfield_default_local(thermal,'hA_os', hA_base);
    hA_o_env = getfield_default_local(thermal,'hA_o_env', hA_base);
    hA_s_env = getfield_default_local(thermal,'hA_s_env', hA_base);

    V_oil_per = resFactor_loc * (Ap_single * (2*Lgap));
    n_tot = sum(multi_col);
    m_oil_tot = n_tot * (rho_base * V_oil_per);
    m_steel_tot = steel_to_oil_mass_ratio * m_oil_tot;
    C_oil = max(m_oil_tot * cp_oil_loc, eps);
    C_steel = max(m_steel_tot * cp_steel_loc, eps);

    z0 = zeros(2*n + 2,1);
    z0(2*n + 1) = T0_C;
    z0(2*n + 2) = T0_C;

    if numel(t) > 1
        dt_nom = median(diff(t));
    else
        dt_nom = 1;
    end
    dt_nom = max(dt_nom, eps);

    opts = odeset('RelTol',1e-3,'AbsTol',1e-6, ...
                  'MaxStep', 0.05);

    ctx = struct('n',n,'Ns',Ns,'M',M,'C',C,'K',K,'r',r,'agf',agf, ...
                 'Nvec',Nvec,'Mvec',Mvec,'multi',multi_row,'multi_col',multi_col,'k_sd',k_sd,'c_lam',c_lam, ...
                 'Ap_story',Ap_story_row,'Ap_story_col',Ap_story_col,'Ao_story',Ao_story_row,'Ao_story_col',Ao_story_col,'Ap_single',Ap_single, ...
                 'Qcap_story',Qcap_story,'Qcap_story_col',Qcap_story.','n_orf_total',n_orf_total_row,'n_orf_total_col',n_orf_total_col,'d_o_single',d_o_single, ...
                 'Lori',Lori,'mu_ref',mu_ref,'b_mu',b_mu,'mu_floor_flag',logical(mu_floor_flag), ...
                 'mu_min_phys',mu_min_phys,'mu_max_phys',mu_max_phys,'gain_pf',gain_pf,'pf_res_only',logical(pf_res_only), ...
                 'cfg',cfg,'soft_eps',soft_eps,'rho_base',rho_base,'alpha_rho',alpha_rho, ...
                 'beta0',beta0,'b_beta',b_beta,'T_ref',T_ref,'T_env',T_env, ...
                 'C_oil',C_oil,'C_steel',C_steel,'hA_os',hA_os,'hA_o_env',hA_o_env,'hA_s_env',hA_s_env, ...
                 'T_min',T0_C,'T_max',T0_C + dT_max,'orf',orf,'Qcap_eff',Qcap_eff, ...
                 'thermal',thermal, 'rho_min', getfield_default_local(thermal,'rho_min',600), 'gate_k', getfield_default_local(cfg,'gate_k',20), ...
                 'hyd_inertia', logical(getfield_default_local(cfg_on,'hyd_inertia', false)), 'dt', dt_nom, ...
                 'cache_token', randi([0 2^31-1]));

    try
        sol = ode15s(@(tt,z) damper_rhs_local(tt,z,ctx), [t(1) t(end)], z0, opts);
    catch
        opts2 = odeset(opts,'RelTol',5e-3,'AbsTol',2e-6);
        sol = ode23tb(@(tt,z) damper_rhs_local(tt,z,ctx), [t(1) t(end)], z0, opts2);
    end

    t_end_num = sol.x(end);
    t_clip = min(t, t_end_num);
    Z = deval(sol, t_clip).';

    x_eff = Z(:,1:n);
    v_eff = Z(:,n+1:2*n);
    T_o = Z(:,2*n+1);
    T_s = Z(:,2*n+2);

    diag_eff = postprocess_state_local(ctx, t_clip, x_eff, v_eff, T_o, T_s);

    Nt_eff = numel(t_clip);
    ag_clip = agf(t_clip);

    F_story = diag_eff.story_force;
    F = zeros(Nt_eff, n);
    for k = 1:Ns
        idxN = Nvec(k);
        idxM = Mvec(k);
        F(:,idxN) = F(:,idxN) - F_story(:,k);
        F(:,idxM) = F(:,idxM) + F_story(:,k);
    end
    a_rel_eff = ( -(M\(C*v_eff.' + K*x_eff.' + F.' )).' - ag_clip(:).*r.' );

    x = x_eff;
    a_rel = a_rel_eff;
    diag = diag_eff;

    if Nt_eff < Nt
        x = pad_timeseries_local(x, Nt_eff, Nt);
        a_rel = pad_timeseries_local(a_rel, Nt_eff, Nt);
        diag = pad_diag_struct_local(diag, Nt_eff, Nt);
    end

    ts = diag;

    function val_pad = pad_timeseries_local(val, Nt_eff_loc, Nt_loc)
        if Nt_eff_loc >= Nt_loc || Nt_loc <= 0
            val_pad = val;
            return;
        end

        pad_len = Nt_loc - Nt_eff_loc;
        if pad_len <= 0 || isempty(val)
            val_pad = val;
            return;
        end

        sz = size(val);
        if sz(1) ~= Nt_eff_loc
            val_pad = val;
            return;
        end

        nd = ndims(val);
        rep_dims = ones(1, max(nd, 1));
        rep_dims(1) = pad_len;

        idx = repmat({':'}, 1, max(nd,1));
        idx{1} = Nt_eff_loc;
        last_slice = val(idx{:});
        pad_block = repmat(last_slice, rep_dims);
        val_pad = cat(1, val, pad_block);
    end

    function out_struct = pad_diag_struct_local(in_struct, Nt_eff_loc, Nt_loc)
        out_struct = in_struct;
        if Nt_eff_loc >= Nt_loc
            return;
        end
        fns = fieldnames(in_struct);
        for jj = 1:numel(fns)
            name = fns{jj};
            val = in_struct.(name);
            if isnumeric(val) || islogical(val)
                sz = size(val);
                if ~isempty(val) && sz(1) == Nt_eff_loc
                    out_struct.(name) = pad_timeseries_local(val, Nt_eff_loc, Nt_loc);
                end
            end
        end
    end

    function dz = damper_rhs_local(tt,z,ctx)
        x_loc = z(1:ctx.n);
        v_loc = z(ctx.n+1:2*ctx.n);
        T_o_loc = z(2*ctx.n+1);
        T_s_loc = z(2*ctx.n+2);

        drift = x_loc(ctx.Mvec) - x_loc(ctx.Nvec);
        dvel = v_loc(ctx.Mvec) - v_loc(ctx.Nvec);

        mu_loc = ctx.mu_ref * exp(ctx.b_mu * (T_o_loc - ctx.T_ref));
        if ctx.mu_floor_flag
            mu_loc = max(ctx.mu_min_phys, mu_loc);
        end
        if isfinite(ctx.mu_max_phys)
            mu_loc = min(ctx.mu_max_phys, mu_loc);
        end
        rho_loc = max(ctx.rho_min, ctx.rho_base / (1 + ctx.alpha_rho * (T_o_loc - ctx.T_ref)));
        beta_loc = max(1e8, ctx.beta0 * exp(ctx.b_beta * (T_o_loc - ctx.T_ref))); %#ok<NASGU>
        p_vap_loc = p_vap_Antoine_local(T_o_loc, ctx.thermal, ctx.orf); %#ok<NASGU>

        F_lin_loc = ctx.k_sd * drift;
        orf_loc = compute_orifice_terms_local(dvel, F_lin_loc, mu_loc, rho_loc, ctx);

        % gate_k zaten ctx içine ekli (default 20); hız ile ölçekli kapı:
        s_loc = tanh(ctx.gate_k * dvel);                     % [-1,1], C^∞ pürüzsüz
        p_pf_loc = - orf_loc.dP_eff .* s_loc;                % dissipatif, süreklilik var

        if ctx.pf_res_only
            p_pf_loc = s_loc .* max(0, s_loc .* p_pf_loc);
        end
        PF_term     = p_pf_loc;
        F_story_loc = F_lin_loc + ctx.Ap_story_col .* PF_term;  % N

        F_loc = zeros(ctx.n,1);
        for ii = 1:ctx.Ns
            F_loc(ctx.Nvec(ii)) = F_loc(ctx.Nvec(ii)) - F_story_loc(ii);
            F_loc(ctx.Mvec(ii)) = F_loc(ctx.Mvec(ii)) + F_story_loc(ii);
        end

        dv_loc = ctx.M \ ( -ctx.C*v_loc - ctx.K*x_loc - F_loc - ctx.M*ctx.r*ctx.agf(tt) );

        qmag_loc = abs(orf_loc.Q);
        P_heat_loc = sum( (orf_loc.dP_eff .* qmag_loc) .* ctx.multi_col );
        P_heat_loc = max(P_heat_loc, 0);

        dT_o_loc = ( P_heat_loc - ctx.hA_os * (T_o_loc - T_s_loc) ...
                     - ctx.hA_o_env * (T_o_loc - ctx.T_env) ) / ctx.C_oil;
        dT_s_loc = ( ctx.hA_os * (T_o_loc - T_s_loc) ...
                     - ctx.hA_s_env * (T_s_loc - ctx.T_env) ) / ctx.C_steel;

        if T_o_loc <= ctx.T_min
            dT_o_loc = max(dT_o_loc, 0);
        elseif T_o_loc >= ctx.T_max
            dT_o_loc = min(dT_o_loc, 0);
        end
        if T_s_loc <= ctx.T_min
            dT_s_loc = max(dT_s_loc, 0);
        elseif T_s_loc >= ctx.T_max
            dT_s_loc = min(dT_s_loc, 0);
        end

        dz = [ v_loc; dv_loc; dT_o_loc; dT_s_loc ];
    end

    function out = postprocess_state_local(ctx, t, x, v, T_o, T_s)
        Nt_loc = numel(t);
        Ns_loc = ctx.Ns;
        drift = x(:,ctx.Mvec) - x(:,ctx.Nvec);
        dvel = v(:,ctx.Mvec) - v(:,ctx.Nvec);

        mu_series = ctx.mu_ref * exp(ctx.b_mu * (T_o - ctx.T_ref));
        if ctx.mu_floor_flag
            mu_series = max(ctx.mu_min_phys, mu_series);
        end
        if isfinite(ctx.mu_max_phys)
            mu_series = min(ctx.mu_max_phys, mu_series);
        end
        rho_series = max(ctx.rho_min, ctx.rho_base ./ (1 + ctx.alpha_rho * (T_o - ctx.T_ref)));
        beta_series = max(1e8, ctx.beta0 * exp(ctx.b_beta * (T_o - ctx.T_ref)));
        p_vap_series = p_vap_Antoine_local(T_o, ctx.thermal, ctx.orf);

        F_lin_series = ctx.k_sd * drift;
        orf_series = compute_orifice_terms_local(dvel, F_lin_series, mu_series, rho_series, ctx);

        % Energy bookkeeping and cavitation mask
        Ap_mat = repmat(ctx.Ap_story, Nt_loc,1);
        s_series = tanh(ctx.gate_k * dvel);
        p_pf_series = - orf_series.dP_eff .* s_series;
        if ctx.pf_res_only
            p_pf_series = s_series .* max(0, s_series .* p_pf_series);
        end
        PF_term_series  = p_pf_series;
        PF_force_series = Ap_mat .* PF_term_series;

        F_story_series = F_lin_series + PF_force_series;

        % --- ORIFICE LOSS: tek kanal ---
        P_loss_series = orf_series.dP_eff .* abs(orf_series.Q);
        P_sum = (P_loss_series) * ctx.multi';

        % --- STRUCTURAL POWER: kapalı (k_sd=0) ---
        P_struct_tot = zeros(Nt_loc,1);

        % --- Enerji ---
        E_orf    = cumtrapz(t, P_loss_series * ctx.multi');
        E_struct = cumtrapz(t, P_struct_tot);

        % Passive safety check (optional logging)
        passivity_viol = any((PF_force_series .* dvel) > 1e-9, 'all');

        p_up_series = ctx.orf.p_amb * ones(size(F_lin_series));
        p_vap_mat = repmat(p_vap_series, 1, Ns_loc);
        p2_est_series = p_up_series - orf_series.dP_eff;
        cav_thresh = ctx.orf.cav_sf * p_vap_mat;
        cav_mask = p2_est_series <= cav_thresh;

        orf_series_diag = struct('dP_eff', orf_series.dP_eff, 'Q', orf_series.Q);

        out = struct('dvel', dvel, 'story_force', F_story_series, 'Q', orf_series.Q, ...
            'dP_eff', orf_series.dP_eff, 'PF', PF_force_series, 'cav_mask', cav_mask, 'P_sum', P_sum, ...
            'E_orf', E_orf, 'E_struct', E_struct, 'passivity_viol', passivity_viol, 'T_oil', T_o, 'T_steel', T_s, ...
            'mu', mu_series, 'rho', rho_series, 'beta', beta_series, 'p_vap', p_vap_series, ...
            'c_lam', ctx.c_lam, 'dP_kv', orf_series.dP_kv, 'dP_cav', orf_series.dP_cav, ...
            'F_lin', F_lin_series, 'Re', orf_series.Re, 'orf_series', orf_series_diag);

    end
    function out = compute_orifice_terms_local(dvel, F_lin, mu_val, rho_val, ctx)
        Ns_loc = ctx.Ns;
        isColumn = isvector(dvel);
        if isColumn
            dvel_mat = reshape(dvel(:), 1, Ns_loc);
            F_lin_mat = reshape(F_lin(:), 1, Ns_loc);
            mu_mat = expand_to_matrix_local(mu_val, 1, Ns_loc);
            rho_mat = expand_to_matrix_local(rho_val, 1, Ns_loc);
        else
            dvel_mat = dvel;
            if size(dvel_mat,2) ~= Ns_loc
                error('compute_orifice_terms_local:bad_dvel','dvel column count mismatch');
            end
            F_lin_mat = F_lin;
            mu_mat = expand_to_matrix_local(mu_val, size(dvel_mat,1), Ns_loc);
            rho_mat = expand_to_matrix_local(rho_val, size(dvel_mat,1), Ns_loc);
        end

        nRows = size(dvel_mat,1);
        Ap_mat    = repmat(ctx.Ap_story,   nRows, 1);
        Ao_mat    = repmat(ctx.Ao_story,   nRows, 1);
        n_orf_mat = repmat(ctx.n_orf_total, nRows, 1);

        persistent cache_tokens cache_lastQ
        if isempty(cache_tokens)
            cache_tokens = zeros(0,1);
            cache_lastQ = {};
        end

        cache_token = NaN;
        if isfield(ctx,'cache_token')
            cache_token = double(ctx.cache_token);
        end

        Q = Ap_mat .* dvel_mat;                   % (1) Kinematic continuity

        dt_loc = getfield_default_local(ctx,'dt',0);
        if size(Q,1) > 1 && dt_loc > 0
            dQdt = gradient(Q, dt_loc);
            if ~isnan(cache_token)
                idx_tok = find(cache_tokens == cache_token, 1);
                if isempty(idx_tok)
                    cache_tokens(end+1,1) = cache_token; %#ok<AGROW>
                    cache_lastQ{end+1} = Q(end,:);
                else
                    cache_lastQ{idx_tok} = Q(end,:);
                end
            end
        elseif dt_loc > 0 && ~isnan(cache_token)
            idx_tok = find(cache_tokens == cache_token, 1);
            if isempty(idx_tok)
                cache_tokens(end+1,1) = cache_token; %#ok<AGROW>
                cache_lastQ{end+1} = Q;
                dQdt = zeros(size(Q));
            else
                prev_Q = cache_lastQ{idx_tok};
                if isempty(prev_Q)
                    dQdt = zeros(size(Q));
                else
                    dQdt = (Q - prev_Q) ./ dt_loc;
                end
                cache_lastQ{idx_tok} = Q;
            end
        else
            dQdt = zeros(size(Q));
        end

        Re = (rho_mat .* abs(Q)) .* ctx.d_o_single ./ max(Ao_mat .* mu_mat, 1e-9);
        Cd = ctx.orf.CdInf - (ctx.orf.CdInf - ctx.orf.Cd0) ./ (1 + Re.^ctx.orf.p_exp);
        Cd_min = getfield_default_local(ctx.orf,'Cd_min',0.5); Cd_max = getfield_default_local(ctx.orf,'Cd_max',0.9); Cd = min(max(Cd, Cd_min), Cd_max);

        dP_kv = 0.5 * rho_mat .* (abs(Q) ./ max(Cd .* Ao_mat, 1e-12)).^2;
        R_lam = (128 * mu_mat .* ctx.Lori) ./ max(pi * ctx.d_o_single^4, 1e-24) ./ max(n_orf_mat,1);
dP_lam = R_lam .* Q;    % ← işaretli laminer düşüm
        dP_h = dP_kv + dP_lam;
        if isfield(ctx,'hyd_inertia') && ctx.hyd_inertia
            Lori_mat = expand_to_matrix_local(ctx.Lori, size(Q,1), size(Q,2));
            Lh1 = rho_mat .* Lori_mat ./ max(Ao_mat, 1e-12);
            Lh = Lh1 ./ max(n_orf_mat, 1);
            dP_h = dP_h + Lh .* abs(dQdt);
        end

        p_up = ctx.orf.p_amb * ones(size(dP_h));
        p_vap = getfield_default_local(ctx.orf,'p_vap_eff', 150);
        p_vap_mat = expand_to_matrix_local(p_vap, size(dP_h,1), size(dP_h,2));
        gamma_c = getfield_default_local(ctx.orf,'cav_sf', 1.0);
        dP_cav = max(p_up - gamma_c * p_vap_mat, 0);

        dP_full = dP_h;
        epsm = max(getfield_default_local(ctx,'soft_eps',1e5), ...
                   0.05*median(dP_full(:) + 1));
        dP_eff = 0.5*(dP_full + dP_cav - sqrt((dP_full - dP_cav).^2 + epsm.^2));

        sgn   = sign(Q + 0);
F_orf = - dP_eff .* Ap_mat .* sgn; % ← hızın tersine (dissipatif)


        if isColumn
            F_orf = F_orf.';
            Q = Q.';
            dP_kv = dP_kv.';
            dP_lam = dP_lam.';
            dP_h = dP_h.';
            dP_cav = dP_cav.';
            dP_eff = dP_eff.';
            Re = Re.';
            Cd = Cd.';
        end

        out = struct('F_orf',F_orf,'Q',Q,'dP_kv',dP_kv,'dP_lam',dP_lam,'dP_h',dP_h, ...
            'dP_cav',dP_cav,'dP_eff',dP_eff,'Re',Re,'Cd',Cd);
    end

end
function mat = expand_to_matrix_local(val, nR, nC)
    if isscalar(val)
        mat = val * ones(nR, nC);
    else
        sz = size(val);
        if isequal(sz, [nR, nC])
            mat = val;
        elseif isvector(val) && numel(val) == nC
            mat = repmat(val(:).', nR, 1);
        elseif isvector(val) && numel(val) == nR
            mat = repmat(val(:), 1, nC);
        else
            error('expand_to_matrix_local:bad_size','Unsupported input dimensions');
        end
    end
end
function p_v = p_vap_Antoine_local(T_C, thermal, orf)
    %#ok<INUSD>
    p_v_const = getfield_default_local(thermal,'p_vap_const_Pa', 150);
    p_v = max(1, p_v_const);
end
function params = default_params()

    n  = 10;
    m  = 360e3 * ones(n,1);
    k  = 6.5e8 * ones(n,1);
    c  = 6.2e6 * ones(n,1);
    story_height = 3.0;

    M  = diag(m);
    K  = zeros(n); C0 = zeros(n);
    for ii = 1:n
        kL = k(ii); cL = c(ii);
        if ii < n
            kU = k(ii+1); cU = c(ii+1);
        else
            kU = 0; cU = 0;
        end
        K(ii,ii) = kL + kU; C0(ii,ii) = cL + cU;
        if ii > 1
            K(ii,ii-1) = -kL; C0(ii,ii-1) = -cL;
        end
        if ii < n
            K(ii,ii+1) = -kU; C0(ii,ii+1) = -cU;
        end
    end

    [~,D] = eig(K,M);
    w  = sqrt(sort(diag(D),'ascend'));
    T1 = 2*pi/w(1);

    params = struct();
    params.n = n;
    params.M = M; params.C0 = C0; params.K = K;
    params.story_height = story_height;
    params.T1 = T1;

    params.Dp   = 0.125;
    params.Lgap = 0.055;
    params.orf = struct('Cd0',0.61,'CdInf',0.80,'Rec',3000,'p_exp',1.1, ...
        'p_amb',1.0e5,'p_cav_eff',2.0e3,'cav_sf',1.0,'d_o',3.0e-3,'veps',0.10);
    params.Lori = 0.10;
    params.mu_ref = 0.9;

    params.Kd = 1.6e9;
    params.Ebody = 2.1e11;
    params.Gsh = 79e9;
    params.d_w = 12e-3;
    params.D_m = 80e-3;
    params.n_turn = 8;

    params.rho = 850;
    params.n_orf = 6;

    params.T0_C = 25;
    params.T_ref_C = 25;
    params.b_mu = -0.013;
    params.thermal = struct('hA_W_perK',450,'hA_os',450,'hA_o_env',450,'hA_s_env',300, ...
        'T_env_C',25,'max_iter',3,'tol_K',0.5,'relax',0.5,'dT_max',80,'Vmin_fac',0.97, ...
        'rho_ref',850,'alpha_rho',7e-4,'beta0',1.6e9,'b_beta',-0.0035, ...
        'antoine_A',5.0,'antoine_B',1700,'antoine_C',-80);


    params.steel_to_oil_mass_ratio = 1.5;
    params.n_dampers_per_story = 2;
    params.story_mask = ones(n-1,1);
    params.cp_oil = 1800;
    params.cp_steel = 500;
    params.resFactor = 12;

    params.c_lam_cap = 2e7;
    params.c_lam_min_frac = 0.05;
    params.c_lam_min_abs = 1e5;

    cfg = struct();
    cfg.on = struct('pressure_force',true,'mu_floor',false, ...
        'pf_resistive_only',false,'Rlam',true,'Rkv',true, ...
        'hyd_inertia',false,'cavitation',true);
    cfg.use_orifice = true;
    cfg.use_thermal = true;
    cfg.K_leak = 0;
    cfg.compat_simple = false;
    cfg.num = struct('softmin_eps',1e5,'mu_min_phys',0.6, ...
        'dP_cap',30e6,'Qcap_scale',1);

    params.cfg = cfg;

    params = build_params_local(params);
end

function params = build_params_local(params)
    params = recompute_damper_params_local(params);
    if isfield(params,'orf') && isfield(params.orf,'CdInf') && ...
            isfield(params,'Ao') && isfield(params,'rho')
        dP_cap = NaN; if isfield(params,'cfg') && isfield(params.cfg,'num'), dP_cap = getfield_default_local(params.cfg.num,'dP_cap',NaN); end; if isfinite(dP_cap)
        params.Qcap_big = max(params.orf.CdInf * params.Ao, 1e-9) * sqrt(2*dP_cap / params.rho);
    else
        params.Qcap_big = 0;
    end;
    else
        params.Qcap_big = 0;
    end
    if isfield(params,'c_lam_min_abs') && isfield(params,'c_lam_min_frac') && ...
            isfield(params,'c_lam0')
        params.c_lam_min = max(params.c_lam_min_abs, ...
            params.c_lam_min_frac * params.c_lam0);
    end
end

function params = recompute_damper_params_local(params)
    if ~isstruct(params), return; end

    if isfield(params,'Dp_mm'),    params.Dp   = params.Dp_mm/1000; end
    if isfield(params,'d_w_mm'),   params.d_w  = params.d_w_mm/1000; end
    if isfield(params,'D_m_mm'),   params.D_m  = params.D_m_mm/1000; end
    if isfield(params,'Lori_mm'),  params.Lori = params.Lori_mm/1000; end
    if isfield(params,'orf') && isfield(params.orf,'d_o_mm')
        params.orf.d_o = params.orf.d_o_mm/1000;
    end

    req = {'Dp','d_w','D_m','n_turn','mu_ref','Lori','Lgap','Kd','Ebody','Gsh'};
    if ~all(isfield(params,req)) || ~isfield(params,'orf') || ~isfield(params.orf,'d_o')
        return;
    end

    if ~isfield(params,'n_orf') && isfield(params,'orf') && isfield(params.orf,'n_orf')
        params.n_orf = params.orf.n_orf;
    end
    if ~isfield(params,'n_orf'), params.n_orf = 1; end

    Ap = pi * params.Dp^2 / 4;
    Ao_single = pi * params.orf.d_o^2 / 4;
    Ao = params.n_orf * Ao_single;

    k_h = params.Kd * Ap^2 / params.Lgap;
    k_s = params.Ebody * Ap / params.Lgap;
    k_hyd = 1 / (1/k_h + 1/k_s);
    k_p = params.Gsh * params.d_w^4 / (8 * params.n_turn * params.D_m^3);

    % --- Laminer eşdeğer c: R_lam * A_p^2 ile TUTARLI ---
    R_lam_single = (128 * params.mu_ref * params.Lori) / (pi * params.orf.d_o^4);
    R_lam_total  = R_lam_single / max(params.n_orf,1);
    c_lam0 = R_lam_total * ((pi * params.Dp^2 / 4)^2);

    params.Ap = Ap;
    params.Ao = Ao;

    % Disable hydraulic and parallel spring stiffness contributions
    params.k_p = 0;
    params.k_sd = 0;
    params.c_lam0 = c_lam0;
end

function val = getfield_default_local(S, fname, defaultVal)
    if ~isstruct(S) || ~isfield(S,fname) || isempty(S.(fname))
        val = defaultVal; return;
    end
    valRaw = S.(fname);
    if isnumeric(valRaw)
        if isempty(valRaw) || any(~isfinite(valRaw(:)))
            val = defaultVal;
        else
            val = valRaw;
        end
    else
        val = valRaw;
    end
end
function metrics = compute_penalties_local(metrics)
    lambda = 5;
    pwr = 1.0;
    W = struct('dP',1,'Qcap',1,'cav',1.5,'T',0.5,'mu',0.5);
    thr = default_qc_thresholds_local(struct());

    metrics.pen_dP = penalty_rel_local(metrics.dP95, thr.dP95_max, pwr);
    metrics.pen_Qcap = penalty_rel_local(metrics.Qcap95, thr.Qcap95_max, pwr);

    cav_val = metrics.cav_pct;
    if ~isfinite(cav_val)
        cav_val = 0;
    end
    if thr.cav_pct_max <= 0
        metrics.pen_cav = max(0, cav_val)^pwr;
    else
        metrics.pen_cav = penalty_rel_local(cav_val, thr.cav_pct_max, pwr);
    end

    metrics.pen_T = penalty_rel_local(metrics.T_end, thr.T_end_max, pwr);
    metrics.pen_mu = penalty_rev_local(metrics.mu_end, thr.mu_end_min, pwr);

    metrics.pen = lambda * (W.dP*metrics.pen_dP + W.Qcap*metrics.pen_Qcap + ...
                            W.cav*metrics.pen_cav + W.T*metrics.pen_T + W.mu*metrics.pen_mu);
end

function val = penalty_rel_local(v, lim, pwr)
    if ~isfinite(v) || ~isfinite(lim)
        val = 0;
        return;
    end
    denom = max(lim, eps);
    val = max(0, (v - lim) / denom)^pwr;
end

function val = penalty_rev_local(v, lim, pwr)
    if ~isfinite(v) || ~isfinite(lim)
        val = 0;
        return;
    end
    denom = max(lim, eps);
    val = max(0, (lim - v) / denom)^pwr;
end

function thr = default_qc_thresholds_local(optsThr)
    if nargin < 1 || isempty(optsThr)
        optsThr = struct();
    end
    thr = struct();
    thr.dP95_max   = getfield_default_local(optsThr,'dP95_max',50e6);
    thr.Qcap95_max = getfield_default_local(optsThr,'Qcap95_max',0.5);
    thr.cav_pct_max= getfield_default_local(optsThr,'cav_pct_max',0);
    thr.T_end_max  = getfield_default_local(optsThr,'T_end_max',75);
    thr.mu_end_min = getfield_default_local(optsThr,'mu_end_min',0.5);
    extra = setdiff(fieldnames(optsThr), fieldnames(thr));
    for ii = 1:numel(extra)
        thr.(extra{ii}) = optsThr.(extra{ii});
    end
end











