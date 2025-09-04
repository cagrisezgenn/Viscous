function summary = run_batch_windowed(scaled, params, opts)
%RUN_BATCH_WINDOWED Analyse multiple records with windowed metrics.
%   SUMMARY = RUN_BATCH_WINDOWED(SCALED, PARAMS, OPTS) processes each
%   groundâmotion record in the struct array SCALED using
%   RUN_ONE_RECORD_WINDOWED and returns a summary table of key metrics.
%   PARAMS bundles structural and damper properties. OPTS are forwarded to
%   RUN_ONE_RECORD_WINDOWED.
%
%   QC logs are printed for IM consistency, low Arias coverage, physical
%   plausibility of response metrics and saturation/cavitation checks.
%   After processing all records, the worst peak floor acceleration and
%   inter-story drift ratio are reported.

if nargin < 3, opts = struct(); end
n = numel(scaled);

% containers for detailed outputs
all_out = cell(n,1);

% preallocate summary arrays
names      = cell(n,1);
scale      = zeros(n,1);
SaT1       = zeros(n,1);
t5         = zeros(n,1);
t95        = zeros(n,1);
coverage   = zeros(n,1);
PFA_top    = zeros(n,1);
IDR_max    = zeros(n,1);
dP95       = zeros(n,1);
Qcap95     = zeros(n,1);
cav_pct    = zeros(n,1);
T_end      = zeros(n,1);
mu_end     = zeros(n,1);
E_orf_win  = zeros(n,1);
E_struct_win = zeros(n,1);
E_ratio_win = zeros(n,1);
zeta1_hot  = zeros(n,1);
z2_over_z1 = zeros(n,1);
pass_fail  = cell(n,1);

worstPFA = -inf; worstPFA_name = '';
worstIDR = -inf; worstIDR_name = '';

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

    m = out.metr;
    PFA_top(k)  = m.PFA_top;
    IDR_max(k)  = m.IDR_max;
    dP95(k)     = m.dP_orf_q95;
    Qcap95(k)   = m.Qcap_ratio_q95;
    cav_pct(k)  = m.cav_pct;
    T_end(k)    = m.T_oil_end;
    mu_end(k)   = m.mu_end;
    E_orf_win(k)   = m.E_orifice_win;
    E_struct_win(k)= m.E_struct_win;
    E_ratio_win(k) = m.E_ratio_win;
    zeta1_hot(k)   = m.zeta1_hot;
    if isfield(m,'z2_over_z1_hot')
        z2_over_z1(k) = m.z2_over_z1_hot;
    else
        z2_over_z1(k) = NaN;
    end

    %% ------------------------- QC logs ---------------------------
    fail = false;
    if ~isfinite(SaT1(k)) || SaT1(k) <= 0
        fprintf('[QC] %s: IM inconsistency (SaT1=%.3f)\n', names{k}, SaT1(k));
        fail = true;
    end
    if out.win.flag_low_arias
        fprintf('[QC] %s: low Arias coverage (%.3f)\n', names{k}, coverage(k));
        fail = true;
    end
    if ~isfinite(PFA_top(k)) || ~isfinite(IDR_max(k))
        fprintf('[QC] %s: physical metrics invalid (PFA=%.3f, IDR=%.3f)\n', ...
            names{k}, PFA_top(k), IDR_max(k));
        fail = true;
    end
    if Qcap95(k) > 1
        fprintf('[QC] %s: saturation (Qcap95=%.3f)\n', names{k}, Qcap95(k));
        fail = true;
    end
    if cav_pct(k) > 0
        fprintf('[QC] %s: cavitation %.1f%%\n', names{k}, 100*cav_pct(k));
        fail = true;
    end
    if fail
        pass_fail{k} = 'fail';
    else
        pass_fail{k} = 'pass';
    end

    if PFA_top(k) > worstPFA
        worstPFA = PFA_top(k);
        worstPFA_name = names{k};
    end
    if IDR_max(k) > worstIDR
        worstIDR = IDR_max(k);
        worstIDR_name = names{k};
    end
end

% assemble summary table
summary = table(names, scale, SaT1, t5, t95, coverage, PFA_top, IDR_max, ...
    dP95, Qcap95, cav_pct, T_end, mu_end, E_orf_win, E_struct_win, ...
    E_ratio_win, zeta1_hot, z2_over_z1, pass_fail, ...
    'VariableNames', {'name','scale','SaT1','t5','t95','coverage', ...
    'PFA_top','IDR_max','dP95','Qcap95','cav_pct','T_end','mu_end', ...
    'E_orf_win','E_struct_win','E_ratio_win','zeta1_hot','z2_over_z1', ...
    'pass_fail'});

fprintf('Worst PFA: %s (%.3f)\n', worstPFA_name, worstPFA);
fprintf('Worst IDR: %s (%.3f)\n', worstIDR_name, worstIDR);

end
