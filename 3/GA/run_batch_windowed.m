function [summary, all_out] = run_batch_windowed(scaled, params, opts)
%RUN_BATCH_WINDOWED Pencereli metriklerle birden fazla kaydı analiz eder.
%   [SUMMARY, ALL_OUT] = RUN_BATCH_WINDOWED(SCALED, PARAMS, OPTS) fonksiyonu,
%   SCALED yapı dizisindeki her yer hareketi kaydını
%   RUN_ONE_RECORD_WINDOWED ile işler ve temel metriklerin özet tablosunu
%   döndürür. ALL_OUT hücre dizisi her kayıt için tam çıktıları içerir.
%   PARAMS, yapısal ve damper özelliklerini; OPTS ise
%   RUN_ONE_RECORD_WINDOWED'e iletilen ayarları barındırır.
%
%   IM tutarlılığı, düşük Arias kapsamı, tepki metriklerinin fiziksel
%   geçerliliği ve satürasyon/kavitasyon kontrolleri için QC logları
%   yazdırılır. Tüm kayıtlar işlendiğinde, en kötü tepe kat ivmesi ve
%   katlar arası ötelenme oranı ile bunlara karşılık gelen \mu faktörü
%   raporlanır.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'mu_factors'), opts.mu_factors = [0.75 1.00 1.25]; end
if ~isfield(opts,'mu_weights'), opts.mu_weights = [0.2 0.6 0.2]; end
if ~isfield(opts,'thr'), opts.thr = struct(); end
opts.thr = Utils.default_qc_thresholds(opts.thr);

assert(isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), ...
    'run_batch_windowed: params.thermal.hA_W_perK eksik');

do_export = isfield(opts,'do_export') && opts.do_export;
if do_export
    if isfield(opts,'outdir')
        outdir = opts.outdir;
    else
        ts = datestr(now,'yyyymmdd_HHMMSS_FFF');
        outdir = fullfile('out', ts);
    end
    if ~exist(outdir,'dir'), mkdir(outdir); end
    if ~Utils.getfield_default(opts,'quiet',false)
        diary(fullfile(outdir,'console.log'));
    end
else
    outdir = '';
end

%% Girdi Hazırlığı
% Çalışma için gerekli dizilerin hazırlanması
n = numel(scaled);

all_out = cell(n,1);

names    = cell(n,1);
scale    = zeros(n,1);
SaT1     = zeros(n,1);
t5       = zeros(n,1);
t95      = zeros(n,1);
coverage = zeros(n,1);
rank_score = nan(n,1);

    % politika/sıra bilgisi
policy_val = Utils.getfield_default(opts,'thermal_reset','each');
order_val  = Utils.getfield_default(opts,'order','natural');
policy_col = repmat({policy_val}, n,1);
order_col  = repmat({order_val}, n,1);
if isfield(opts,'cooldown_s')
    cooldown_val = opts.cooldown_s;
else
    cooldown_val = NaN;
end
cooldown_col = repmat(cooldown_val, n,1);

PFA_nom    = zeros(n,1);
IDR_nom    = zeros(n,1);
dP95_nom   = zeros(n,1);
Qcap95_nom = zeros(n,1);
cav_nom    = zeros(n,1);
zeta1_hot       = zeros(n,1);
z2_over_z1_hot  = zeros(n,1);
P_mech      = zeros(n,1);
Re_max      = zeros(n,1);

PFA_w    = zeros(n,1);
IDR_w    = zeros(n,1);
dP95_w   = zeros(n,1);
Qcap95_w = zeros(n,1);
Q_q95_w  = zeros(n,1);
Q_q50_w  = zeros(n,1);
dP50_w   = zeros(n,1);

PFA_worst    = zeros(n,1);
IDR_worst    = zeros(n,1);
dP95_worst   = zeros(n,1);
Qcap95_worst = zeros(n,1);
Q_q95_worst  = zeros(n,1);
Q_q50_worst  = zeros(n,1);
dP_orf_q50_worst = zeros(n,1);
PF_p95_worst = zeros(n,1);
which_mu_PFA = zeros(n,1);
which_mu_IDR = zeros(n,1);
T_end_worst  = zeros(n,1);
mu_end_worst = zeros(n,1);
cav_pct_worst = zeros(n,1);
    % Damperli üst kat için yeni en kötü durum tepe metrikleri
x10_max_D_worst = zeros(n,1);
a10abs_max_D_worst = zeros(n,1);
    % Enerji özetleri (mu boyunca en kötü değer)
E_orifice_sum = zeros(n,1);
E_struct_sum  = zeros(n,1);
E_ratio       = zeros(n,1);
qc_all_mu    = false(n,1);

T_start    = zeros(n,1);
T_end      = zeros(n,1);
mu_end     = zeros(n,1);
clamp_hits = zeros(n,1);

worstPFA = -inf; worstPFA_name = ''; worstPFA_mu = NaN;
worstIDR = -inf; worstIDR_name = ''; worstIDR_mu = NaN;

%% Kayıt Döngüsü
% Her kayıt için pencere analizi
prev_diag = [];
for k = 1:n
    rec = scaled(k);
    out = run_one_record_windowed(rec, params, opts, prev_diag);
    prev_diag = out.diag;
    all_out{k} = out; %#ok<AGROW>

    names{k}    = out.name;
    scale(k)    = out.scale;
    SaT1(k)     = out.SaT1;
    t5(k)       = out.win.t5;
    t95(k)      = out.win.t95;
    coverage(k) = out.win.coverage;
    % rank skoru yalnızca order='worst_first' için hesaplanır
    if strcmpi(order_val,'worst_first')
        if isfield(out,'metr') && isfield(out.metr,'E_orifice_win')
            rank_score(k) = out.metr.E_orifice_win;
        end
    end

    T_start(k)    = out.T_start;
    T_end(k)      = out.T_end;
    mu_end(k)     = out.mu_end;
    clamp_hits(k) = out.clamp_hits;

    m_nom = out.metr;
    PFA_nom(k)    = m_nom.PFA_top;
    IDR_nom(k)    = m_nom.IDR_max;
    dP95_nom(k)   = m_nom.dP_orf_q95;
    Qcap95_nom(k) = m_nom.Qcap_ratio_q95;
    cav_nom(k)    = m_nom.cav_pct;
    zeta1_hot(k)       = Utils.getfield_default(m_nom,'zeta1_hot',NaN);
    z2_over_z1_hot(k)  = Utils.getfield_default(m_nom,'z2_over_z1_hot',NaN);
    P_mech(k)          = Utils.getfield_default(m_nom,'P_mech',NaN);
    Re_max(k)          = Utils.getfield_default(m_nom,'Re_max',NaN);

    m_w = out.weighted;
    PFA_w(k)    = m_w.PFA_top;
    IDR_w(k)    = m_w.IDR_max;
    dP95_w(k)   = m_w.dP_orf_q95;
    Qcap95_w(k) = m_w.Qcap_ratio_q95;
    fields  = {'Q_q95','Q_q50','dP_orf_q50'};
    targets = {Q_q95_w, Q_q50_w, dP50_w};
    for f = 1:numel(fields)
        targets{f} = Utils.assign_if_field(m_w, fields{f}, targets{f}, k);
    end
    [Q_q95_w, Q_q50_w, dP50_w] = targets{:};

    m_ws = out.worst;
    PFA_worst(k)    = m_ws.PFA_top;
    IDR_worst(k)    = m_ws.IDR_max;
    dP95_worst(k)   = m_ws.dP_orf_q95;
    Qcap95_worst(k) = m_ws.Qcap_ratio_q95;
    fields  = {'Q_q95','Q_q50','dP_orf_q50','PF_p95'};
    targets = {Q_q95_worst, Q_q50_worst, dP_orf_q50_worst, PF_p95_worst};
    for f = 1:numel(fields)
        targets{f} = Utils.assign_if_field(m_ws, fields{f}, targets{f}, k);
    end
    [Q_q95_worst, Q_q50_worst, dP_orf_q50_worst, PF_p95_worst] = targets{:};
    which_mu_PFA(k) = m_ws.which_mu.PFA_top;
    which_mu_IDR(k) = m_ws.which_mu.IDR_max;
    T_end_worst(k)  = m_ws.T_oil_end;
    mu_end_worst(k) = m_ws.mu_end;
    cav_pct_worst(k)= m_ws.cav_pct;
    % Üst kat zarf metrikleri
    fields  = {'x10_max_D','a10abs_max_D'};
    targets = {x10_max_D_worst, a10abs_max_D_worst};
    for f = 1:numel(fields)
        targets{f} = Utils.assign_if_field(m_ws, fields{f}, targets{f}, k);
        if ~isfield(m_ws, fields{f})
            fallback = strrep(fields{f}, '_max_', '_pk_');
            targets{f} = Utils.assign_if_field(m_ws, fallback, targets{f}, k);
        end
    end
    [x10_max_D_worst, a10abs_max_D_worst] = targets{:};
    % Enerji toplamları
    fields  = {'E_orifice_full','E_struct_full','E_ratio_full'};
    targets = {E_orifice_sum, E_struct_sum, E_ratio};
    for f = 1:numel(fields)
        targets{f} = Utils.assign_if_field(m_ws, fields{f}, targets{f}, k);
    end
    [E_orifice_sum, E_struct_sum, E_ratio] = targets{:};
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

%% Özet Tablo
% Hesaplanan metrikleri tabloya dönüştür
summary = struct();

summary.table = table(names, scale, SaT1, t5, t95, coverage, rank_score, policy_col, order_col, cooldown_col, ...
    PFA_nom, IDR_nom, dP95_nom, Qcap95_nom, cav_nom, zeta1_hot, z2_over_z1_hot, P_mech, Re_max, ...
    PFA_w, IDR_w, dP95_w, Qcap95_w, Q_q95_w, Q_q50_w, dP50_w, ...
    PFA_worst, IDR_worst, dP95_worst, dP_orf_q50_worst, Qcap95_worst, Q_q95_worst, Q_q50_worst, PF_p95_worst, ...
    cav_pct_worst, x10_max_D_worst, a10abs_max_D_worst, E_orifice_sum, E_struct_sum, E_ratio, ...
    which_mu_PFA, which_mu_IDR, T_end_worst, mu_end_worst, qc_all_mu, ...
    T_start, T_end, mu_end, clamp_hits, ...
    'VariableNames', {'name','scale','SaT1','t5','t95','coverage','rank_score','policy','order','cooldown_s', ...
    'PFA_nom','IDR_nom','dP95_nom','Qcap95_nom','cav_nom','zeta1_hot','z2_over_z1_hot','P_mech','Re_max', ...
    'PFA_w','IDR_w','dP95_w','Qcap95_w','Q_q95_w','Q_q50_w','dP50_w', ...
    'PFA_worst','IDR_worst','dP95_worst','dP_orf_q50_worst','Qcap95_worst','Q_q95_worst','Q_q50_worst','PF_p95_worst', ...
    'cav_pct_worst','x10_max_D_worst','a10abs_max_D_worst','E_orifice_sum','E_struct_sum','E_ratio', ...
    'which_mu_PFA','which_mu_IDR','T_end_worst','mu_end_worst','qc_all_mu', ...
    'T_start','T_end','mu_end','clamp_hits'});

% Tüketicilerle uyumluluk için takma adlar
summary.table.T_oil_end_worst = summary.table.T_end_worst;
summary.table.dP_orf_q95_worst = summary.table.dP95_worst;
try
    summary.table.energy_tot_sum = summary.table.E_orifice_sum + summary.table.E_struct_sum;
catch
end
summary.all_out = all_out;
%% QC Kontrolü
% QC eşiklerine göre sonuçların değerlendirilmesi
% --- summary.csv kullanıcıları için QC bayrakları ve sebep kodları ---
thr = opts.thr;
ok_T    = summary.table.T_end_worst   <= thr.T_end_max;
ok_mu   = summary.table.mu_end_worst  >= thr.mu_end_min;
ok_dP   = summary.table.dP95_worst    <= thr.dP95_max;
ok_Qcap = summary.table.Qcap95_worst  <  thr.Qcap95_max;
ok_cav  = summary.table.cav_pct_worst == 0;
qc_reason = strings(height(summary.table),1);
for r = 1:height(summary.table)
    bad = {};
    if ~ok_T(r),    bad{end+1}='T';  end %#ok<AGROW>
    if ~ok_mu(r),   bad{end+1}='mu'; end %#ok<AGROW>
    if ~ok_dP(r),   bad{end+1}='dP'; end %#ok<AGROW>
    if ~ok_Qcap(r), bad{end+1}='Qcap'; end %#ok<AGROW>
    if ~ok_cav(r),  bad{end+1}='cav'; end %#ok<AGROW>
    qc_reason(r) = strjoin(bad,',');
end
summary.table.ok_T = ok_T;
summary.table.ok_mu = ok_mu;
summary.table.ok_dP = ok_dP;
summary.table.ok_Qcap = ok_Qcap;
summary.table.ok_cav = ok_cav;
summary.table.qc_reason = qc_reason;

if ~isfield(opts,'quiet') || ~opts.quiet
    fprintf('Worst PFA: %s, mu=%.2f\n', worstPFA_name, worstPFA_mu);
    fprintf('Worst IDR: %s, mu=%.2f\n', worstIDR_name, worstIDR_mu);
end

if do_export
    export_results(outdir, scaled, params, opts, summary, all_out);
    if ~Utils.getfield_default(opts,'quiet',false)
        diary off;
    end
end

end
