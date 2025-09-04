function [summary, all_out] = run_batch_windowed(scaled, params, opts)
%RUN_BATCH_WINDOWED Analyse multiple records with windowed metrics.
%   [SUMMARY, ALL_OUT] = RUN_BATCH_WINDOWED(SCALED, PARAMS, OPTS) processes
%   each ground-motion record in the struct array SCALED using
%   RUN_ONE_RECORD_WINDOWED and returns a summary table of key metrics.
%   ALL_OUT is a cell array containing the full output for each record.
%   PARAMS bundles structural and damper properties. OPTS are forwarded to
%   RUN_ONE_RECORD_WINDOWED.
%
%   QC logs are printed for IM consistency, low Arias coverage, physical
%   plausibility of response metrics and saturation/cavitation checks.
%   After processing all records, the worst peak floor acceleration and
%   inter-story drift ratio along with the associated \mu factor are
%   reported.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'mu_factors'), opts.mu_factors = [0.75 1.00 1.25]; end
if ~isfield(opts,'mu_weights'), opts.mu_weights = [0.2 0.6 0.2]; end
n = numel(scaled);

all_out = cell(n,1);

names    = cell(n,1);
scale    = zeros(n,1);
SaT1     = zeros(n,1);
t5       = zeros(n,1);
t95      = zeros(n,1);
coverage = zeros(n,1);

PFA_nom    = zeros(n,1);
IDR_nom    = zeros(n,1);
dP95_nom   = zeros(n,1);
Qcap95_nom = zeros(n,1);
cav_nom    = zeros(n,1);

PFA_w    = zeros(n,1);
IDR_w    = zeros(n,1);
dP95_w   = zeros(n,1);
Qcap95_w = zeros(n,1);

PFA_worst    = zeros(n,1);
IDR_worst    = zeros(n,1);
dP95_worst   = zeros(n,1);
Qcap95_worst = zeros(n,1);
which_mu_PFA = zeros(n,1);
which_mu_IDR = zeros(n,1);
T_end_worst  = zeros(n,1);
mu_end_worst = zeros(n,1);
qc_all_mu    = false(n,1);

worstPFA = -inf; worstPFA_name = ''; worstPFA_mu = NaN;
worstIDR = -inf; worstIDR_name = ''; worstIDR_mu = NaN;

prev_diag = [];
for k = 1:n
    rec = scaled(k);
    out = run_one_record_windowed(rec, [], params, opts, prev_diag);
    prev_diag = out.diag;
    all_out{k} = out; %#ok<AGROW>

    names{k}    = out.name;
    scale(k)    = out.scale;
    SaT1(k)     = out.SaT1;
    t5(k)       = out.win.t5;
    t95(k)      = out.win.t95;
    coverage(k) = out.win.coverage;

    m_nom = out.metr;
    PFA_nom(k)    = m_nom.PFA_top;
    IDR_nom(k)    = m_nom.IDR_max;
    dP95_nom(k)   = m_nom.dP_orf_q95;
    Qcap95_nom(k) = m_nom.Qcap_ratio_q95;
    cav_nom(k)    = m_nom.cav_pct;

    m_w = out.weighted;
    PFA_w(k)    = m_w.PFA_top;
    IDR_w(k)    = m_w.IDR_max;
    dP95_w(k)   = m_w.dP_orf_q95;
    Qcap95_w(k) = m_w.Qcap_ratio_q95;

    m_ws = out.worst;
    PFA_worst(k)    = m_ws.PFA_top;
    IDR_worst(k)    = m_ws.IDR_max;
    dP95_worst(k)   = m_ws.dP_orf_q95;
    Qcap95_worst(k) = m_ws.Qcap_ratio_q95;
    which_mu_PFA(k) = m_ws.which_mu.PFA_top;
    which_mu_IDR(k) = m_ws.which_mu.IDR_max;
    T_end_worst(k)  = m_ws.T_oil_end;
    mu_end_worst(k) = m_ws.mu_end;
    qc_all_mu(k)    = out.qc_all_mu;

    valsPFA = arrayfun(@(s) s.metr.PFA_top, out.mu_results);
    [maxPFA_rec, idxPFA] = max(valsPFA);
    if maxPFA_rec > worstPFA
        worstPFA = maxPFA_rec;
        worstPFA_name = out.name;
        worstPFA_mu   = out.mu_results(idxPFA).mu_factor;
    end
    valsIDR = arrayfun(@(s) s.metr.IDR_max, out.mu_results);
    [maxIDR_rec, idxIDR] = max(valsIDR);
    if maxIDR_rec > worstIDR
        worstIDR = maxIDR_rec;
        worstIDR_name = out.name;
        worstIDR_mu   = out.mu_results(idxIDR).mu_factor;
    end
end

summary = struct();
summary.table = table(names, scale, SaT1, t5, t95, coverage, ...
    PFA_nom, IDR_nom, dP95_nom, Qcap95_nom, cav_nom, ...
    PFA_w, IDR_w, dP95_w, Qcap95_w, ...
    PFA_worst, IDR_worst, dP95_worst, Qcap95_worst, ...
    which_mu_PFA, which_mu_IDR, T_end_worst, mu_end_worst, qc_all_mu, ...
    'VariableNames', {'name','scale','SaT1','t5','t95','coverage', ...
    'PFA_nom','IDR_nom','dP95_nom','Qcap95_nom','cav_nom', ...
    'PFA_w','IDR_w','dP95_w','Qcap95_w', ...
    'PFA_worst','IDR_worst','dP95_worst','Qcap95_worst', ...
    'which_mu_PFA','which_mu_IDR','T_end_worst','mu_end_worst','qc_all_mu'});
summary.all_out = all_out;

fprintf('Worst PFA: %s, mu=%.2f\n', worstPFA_name, worstPFA_mu);
fprintf('Worst IDR: %s, mu=%.2f\n', worstIDR_name, worstIDR_mu);

end
