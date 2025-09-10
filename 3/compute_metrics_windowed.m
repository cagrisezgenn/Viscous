function metr = compute_metrics_windowed(t, x, a_rel, ag, ts, story_height, win, params)
%COMPUTE_METRICS_WINDOWED Compute response metrics in a time window.
%   METR = COMPUTE_METRICS_WINDOWED(T,X,A_REL,AG,TS,STORY_HEIGHT,WIN,PARAMS)
%   evaluates a set of performance metrics for the structural response
%   confined to the time window specified by WIN.IDX.  The fields of METR
%   include peak floor acceleration at the top story, interstory drift
%   ratios, orifice pressure statistics, energy measures and modal damping
%   ratios based on the final ("hot") damper coefficient.
%
%   T, X, A_REL and AG are the time vector, story displacements, relative
%   floor accelerations and ground acceleration.  TS contains additional
%   time-series diagnostics produced by MCK_WITH_DAMPER_TS.  STORY_HEIGHT is
%   the interstory height.  WIN.IDX is a logical vector selecting the window
%   of interest.  PARAMS carries structural and damper parameters including
%   a DIAG field with thermal quantities.

idx = win.idx;

% ---------------- Basic structural response metrics -----------------
% Peak floor acceleration at top story (absolute)
a_top_abs = a_rel(:,end) + ag(:);
metr.PFA_top = max(abs(a_top_abs(idx)));

% === Damperli 10. kat tepe değerleri (tek kayıt/pencere) ===
try
    % x: damperli bağıl yerdeğiştirme [Nt x n], a_top_abs: damperli mutlak ivme [Nt x 1]
    metr.x10_pk_D      = max(abs(x(idx,end)));
    metr.a10abs_pk_D   = max(abs(a_top_abs(idx)));
catch
    if ~isfield(metr,'x10_pk_D'),    metr.x10_pk_D = 0;      end
    if ~isfield(metr,'a10abs_pk_D'), metr.a10abs_pk_D = 0;   end
end

% Top-story absolute displacement (damperli) within window
metr.x10_max_D = max(abs(x(idx,end)));

% Top-story absolute acceleration (damperli) within window
metr.a10abs_max_D = max(abs(a_top_abs(idx)));

% Maximum inter-story drift ratio
drift = (x(:,2:end) - x(:,1:end-1)) / story_height;
metr.IDR_max = max(max(abs(drift(idx,:))));

% --------------- Per-story statistics within the window --------------
q95 = @(A) quantile(A,0.95);
q50 = @(A) quantile(A,0.50);

abs_dP = abs(ts.dP_orf(idx,:));
abs_Q  = abs(ts.Q(idx,:));
Qcap_ratio = ts.Qcap_ratio(idx,:);
abs_story_force = abs(ts.story_force(idx,:));

% 95th percentiles per story
% 50th and 95th percentiles per story
dP_q50         = q50(abs_dP);
dP_q95         = q95(abs_dP);
Q_q50          = q50(abs_Q);
Q_q95          = q95(abs_Q);
Qcap_ratio_q95 = q95(Qcap_ratio);
story_force_q95= q95(abs_story_force);

% Mean cavitation fraction per story
cav_mean = mean(ts.cav_mask(idx,:),1);

% Critical story index based on story force
[metr.story_force_q95, metr.which_story] = max(story_force_q95);
ws = metr.which_story;

% Store corresponding per-story metrics
metr.dP_orf_q95      = dP_q95(ws);
metr.dP_orf_q50      = dP_q50(ws);
metr.Q_q95           = Q_q95(ws);
metr.Q_q50           = Q_q50(ws);
metr.Qcap_ratio_q95  = Qcap_ratio_q95(ws);
metr.cav_pct         = cav_mean(ws);

% ----------------------- Energy calculations -------------------------
w_first = find(idx,1,'first');
w_last  = find(idx,1,'last');
i0 = max(w_first-1,1);

metr.E_orifice_full = ts.E_orf(end);
metr.E_struct_full  = ts.E_struct(end);
metr.E_ratio_full   = metr.E_orifice_full / max(metr.E_struct_full, eps);
% Convenience totals/powers for aggregations
metr.energy_tot = metr.E_orifice_full + metr.E_struct_full;
try
    if isfield(ts,'P_sum') && ~isempty(ts.P_sum)
        metr.P_mech = mean(ts.P_sum(idx));
    end
catch
end

metr.E_orifice_win = ts.E_orf(w_last) - ts.E_orf(i0);
metr.E_struct_win  = ts.E_struct(w_last) - ts.E_struct(i0);
metr.E_ratio_win   = metr.E_orifice_win / max(metr.E_struct_win, eps);

% -------------------- Thermal/viscosity metrics ----------------------
if isfield(params,'diag') && isfield(params.diag,'T_oil')
    metr.T_oil_end = params.diag.T_oil(w_last);
else
    metr.T_oil_end = NaN;
end
if isfield(params,'diag') && isfield(params.diag,'mu')
    metr.mu_end = params.diag.mu(w_last);
else
    metr.mu_end = NaN;
end

% ---------------- Modal damping with hot viscosity -------------------
req_fields = {'M','K','C0','k_sd','toggle_gain','story_mask','n_dampers_per_story'};
if all(isfield(params,req_fields)) && isfield(params,'diag') && isfield(params.diag,'c_lam')
    M  = params.M;  K = params.K;  C0 = params.C0;
    k_sd = params.k_sd;
    Rvec = params.toggle_gain(:);
    nStories = size(M,1) - 1;
    if numel(Rvec)==1, Rvec = Rvec*ones(nStories,1); end
    mask = params.story_mask(:); if numel(mask)==1, mask = mask*ones(nStories,1); end
    ndps = params.n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi = (mask .* ndps);
    c_lam = params.diag.c_lam;

    Kadd = zeros(size(M));
    Cadd = zeros(size(M));
    for i=1:nStories
        idx2 = [i i+1];
        % Use linear R scaling (R), not squared (R^2)
        k_eq = k_sd * Rvec(i) * multi(i);
        c_eq = c_lam * Rvec(i) * multi(i);
        kM = k_eq * [1 -1; -1 1];
        cM = c_eq * [1 -1; -1 1];
        Kadd(idx2,idx2) = Kadd(idx2,idx2) + kM;
        Cadd(idx2,idx2) = Cadd(idx2,idx2) + cM;
    end
    Ktot = K + Kadd;
    Ctot = C0 + Cadd;
    [V,D] = eig(Ktot,M);
    [w2,ord] = sort(diag(D),'ascend');
    V = V(:,ord);
    phi1 = V(:,1); phi2 = V(:,2);
    w1 = sqrt(w2(1)); w2s = sqrt(w2(2));
    n1 = phi1.'*M*phi1; n2 = phi2.'*M*phi2;
    z1 = (phi1.'*Ctot*phi1)/(2*w1*n1);
    z2 = (phi2.'*Ctot*phi2)/(2*w2s*n2);
    metr.zeta1_hot = z1;
    metr.z2_over_z1_hot = z2 / max(z1, eps);
else
    metr.zeta1_hot = NaN;
    metr.z2_over_z1_hot = NaN;
end

% ----------------------- PF metrics (optional) ------------------------
try
    if isfield(ts,'PF')
        PF_abs = abs(ts.PF(idx,:));
        metr.PF_p95 = q95(PF_abs(:,ws));
    end
catch
end

end
