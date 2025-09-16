function [summary, all_out] = run_batch_windowed(scaled, params, opts)
%RUN_BATCH_WINDOWED Pencereli metriklerle birden fazla kaydı analiz eder.
%   [SUMMARY, ALL_OUT] = RUN_BATCH_WINDOWED(SCALED, PARAMS, OPTS) fonksiyonu,
%   SCALED yapı dizisindeki her yer hareketi kaydını
%   RUN_ONE_RECORD_WINDOWED ile işler ve temel metriklerin özet tablosunu
%   döndürür. ALL_OUT hücre dizisi her kayıt için tam çıktıları içerir.
%   PARAMS, yapısal ve damper özelliklerini; OPTS ise
%   RUN_ONE_RECORD_WINDOWED'e iletilen ayarları barındırır.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'thr'), opts.thr = struct(); end
opts.thr = Utils.default_qc_thresholds(opts.thr);

% İstenirse kayıtların yürütülme sırası yeniden düzenlenir
if isfield(opts,'order_perm')
    scaled = scaled(opts.order_perm);
end

assert(isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), ...
    'run_batch_windowed: params.thermal.hA_W_perK eksik');

%% Girdi Hazırlığı
n = numel(scaled);
vars = prepare_inputs(n, params, opts);

%% Kayıt Döngüsü
vars = record_loop(scaled, params, opts, vars);

%% Özet Tablo
summary = build_summary_table(vars, opts);
all_out = vars.all_out;

end

function vars = prepare_inputs(n, params, opts)
%PREPARE_INPUTS Çalışma için gerekli dizileri hazırla
vars = struct();
vars.all_out = cell(n,1);
vars.names    = cell(n,1);
vars.scale    = zeros(n,1);
vars.SaT1     = zeros(n,1);
vars.t5       = zeros(n,1);
vars.t95      = zeros(n,1);
vars.coverage = zeros(n,1);

policy_val = Utils.getfield_default(opts,'thermal_reset','each');
order_val  = Utils.getfield_default(opts,'order','natural');
vars.policy_col = repmat({policy_val}, n,1);
vars.order_col  = repmat({order_val}, n,1);
if isfield(opts,'cooldown_s')
    cooldown_val = opts.cooldown_s;
else
    cooldown_val = NaN;
end
vars.cooldown_col = repmat(cooldown_val, n,1);

vars.PFA    = zeros(n,1);
vars.IDR    = zeros(n,1);
vars.dP95   = zeros(n,1);
vars.Qcap95 = zeros(n,1);
vars.cav_pct = zeros(n,1);
vars.zeta1_hot       = zeros(n,1);
vars.z2_over_z1_hot  = zeros(n,1);
vars.P_mech_sum      = zeros(n,1);
vars.Re_max      = zeros(n,1);
    vars.Q_q95  = zeros(n,1);
    vars.Q_q50  = zeros(n,1);
    vars.dP50   = zeros(n,1);
    vars.PF_p95 = zeros(n,1);
vars.x10_max_damperli = zeros(n,1);
vars.a10abs_max_damperli = zeros(n,1);
vars.E_orifice_sum = zeros(n,1);
vars.E_struct_sum  = zeros(n,1);
vars.energy_tot_sum = zeros(n,1);
vars.E_ratio   = zeros(n,1);
vars.qc_pass   = false(n,1);

vars.T_start    = zeros(n,1);
vars.T_end      = zeros(n,1);
vars.mu_end     = zeros(n,1);
vars.clamp_hits = zeros(n,1);

vars.Dp_mm_col   = repmat(Utils.getfield_default(params,'Dp_mm',NaN), n,1);
vars.mu_ref_col  = repmat(Utils.getfield_default(params,'mu_ref',NaN), n,1);

vars.worstPFA = -inf; vars.worstPFA_name = '';
vars.worstIDR = -inf; vars.worstIDR_name = '';
end

function vars = record_loop(scaled, params, opts, vars)
%RECORD_LOOP Her kaydı pencere analizine tabi tutar
prev_ts = [];
for k = 1:numel(scaled)
    rec = scaled(k);
    out = run_one_record_windowed(rec, params, opts, prev_ts);
    prev_ts = out.ts;
    vars.all_out{k} = out; %#ok<AGROW>

    vars.names{k}    = out.name;
    vars.scale(k)    = out.scale;
    vars.SaT1(k)     = out.SaT1;
    vars.t5(k)       = out.win.t5;
    vars.t95(k)      = out.win.t95;
    vars.coverage(k) = out.win.coverage;

    vars.T_start(k)    = out.T_start;
    vars.T_end(k)      = out.T_end;
    vars.mu_end(k)     = out.mu_end;
    vars.clamp_hits(k) = out.clamp_hits;
    m_nom = out.metr;
    vars.PFA(k)    = m_nom.PFA;
    vars.IDR(k)    = m_nom.IDR;
    vars.dP95(k)   = m_nom.dP95;
    vars.Qcap95(k) = m_nom.Qcap95;
    vars.cav_pct(k)= m_nom.cav_pct;
    vars.zeta1_hot(k)       = Utils.getfield_default(m_nom,'zeta1_hot',NaN);
    vars.z2_over_z1_hot(k)  = Utils.getfield_default(m_nom,'z2_over_z1_hot',NaN);
    vars.P_mech_sum(k)      = Utils.getfield_default(m_nom,'P_mech_sum',NaN);
    vars.Re_max(k)          = Utils.getfield_default(m_nom,'Re_max',NaN);
    vars.Q_q95(k)  = Utils.getfield_default(m_nom,'Q_q95',NaN);
    vars.Q_q50(k)  = Utils.getfield_default(m_nom,'Q_q50',NaN);
    vars.dP50(k)   = Utils.getfield_default(m_nom,'dP50',NaN);
    vars.PF_p95(k) = Utils.getfield_default(m_nom,'PF_p95',NaN);
    vars.x10_max_damperli(k) = Utils.getfield_default(m_nom,'x10_max_damperli',NaN);
    vars.a10abs_max_damperli(k) = Utils.getfield_default(m_nom,'a10abs_max_damperli',NaN);
    vars.E_orifice_sum(k) = Utils.getfield_default(m_nom,'E_orifice_sum',NaN);
    vars.E_struct_sum(k)  = Utils.getfield_default(m_nom,'E_struct_sum',NaN);
    vars.energy_tot_sum(k) = Utils.getfield_default(m_nom,'energy_tot_sum',NaN);
    vars.E_ratio(k)   = Utils.getfield_default(m_nom,'E_ratio',NaN);
    vars.qc_pass(k)   = out.qc_pass;

    if vars.PFA(k) > vars.worstPFA
        vars.worstPFA = vars.PFA(k);
        vars.worstPFA_name = out.name;
    end
    if vars.IDR(k) > vars.worstIDR
        vars.worstIDR = vars.IDR(k);
        vars.worstIDR_name = out.name;
    end
end
end

function summary = build_summary_table(vars, opts)
%BUILD_SUMMARY_TABLE Hesaplanan metrikleri tabloya dönüştür ve QC uygula
summary = struct();

thr = opts.thr;
ok_T    = vars.T_end   <= thr.T_end_max;
ok_mu   = vars.mu_end  >= thr.mu_end_min;
ok_dP   = vars.dP95    <= thr.dP95_max;
ok_Qcap = vars.Qcap95  <  thr.Qcap95_max;
ok_cav  = vars.cav_pct == 0;
qc_reason = strings(numel(vars.names),1);
for r = 1:numel(vars.names)
    bad = {};
    if ~ok_T(r),    bad{end+1}='T';  end %#ok<AGROW>
    if ~ok_mu(r),   bad{end+1}='mu'; end %#ok<AGROW>
    if ~ok_dP(r),   bad{end+1}='dP'; end %#ok<AGROW>
    if ~ok_Qcap(r), bad{end+1}='Qcap'; end %#ok<AGROW>
    if ~ok_cav(r),  bad{end+1}='cav'; end %#ok<AGROW>
    qc_reason(r) = strjoin(bad,',');
end

summary.table = table(vars.names, vars.scale, vars.SaT1, vars.t5, vars.t95, vars.coverage, vars.policy_col, vars.order_col, vars.cooldown_col, ...
    vars.PFA, vars.IDR, vars.dP95, vars.Qcap95, vars.cav_pct, vars.PF_p95, vars.zeta1_hot, vars.z2_over_z1_hot, vars.P_mech_sum, vars.Re_max, ...
    vars.Q_q95, vars.Q_q50, vars.dP50, vars.x10_max_damperli, vars.a10abs_max_damperli, vars.E_orifice_sum, vars.E_struct_sum, vars.energy_tot_sum, vars.E_ratio, vars.qc_pass, ...
    vars.T_start, vars.T_end, vars.mu_end, vars.clamp_hits, vars.Dp_mm_col, vars.mu_ref_col, ...
    ok_T, ok_mu, ok_dP, ok_Qcap, ok_cav, qc_reason, ...
    'VariableNames', {'name','scale','SaT1','t5','t95','coverage','policy','order','cooldown_s', ...
    'PFA','IDR','dP95','Qcap95','cav_pct','PF_p95','zeta1_hot','z2_over_z1_hot','P_mech_sum','Re_max', ...
    'Q_q95','Q_q50','dP50','x10_max_damperli','a10abs_max_damperli','E_orifice_sum','E_struct_sum','energy_tot_sum','E_ratio','qc_pass', ...
    'T_start','T_end','mu_end','clamp_hits','Dp_mm','mu_ref','ok_T','ok_mu','ok_dP','ok_Qcap','ok_cav','qc_reason'});

summary.all_out = vars.all_out;
end

