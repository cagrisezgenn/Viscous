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
if ~isfield(opts,'thr'), opts.thr = struct(); end
opts.thr = Utils.default_qc_thresholds(opts.thr);

% Türetilmiş damper sabitlerini güncelle
params = Utils.recompute_damper_params(params);

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

PFA    = zeros(n,1);
IDR    = zeros(n,1);
dP95   = zeros(n,1);
Qcap95 = zeros(n,1);
cav_pct = zeros(n,1);
zeta1_hot       = zeros(n,1);
z2_over_z1_hot  = zeros(n,1);
P_mech      = zeros(n,1);
Re_max      = zeros(n,1);
Q_q95  = zeros(n,1);
Q_q50  = zeros(n,1);
dP50   = zeros(n,1);
x10_max_D = zeros(n,1);
a10abs_max_D = zeros(n,1);
E_orifice = zeros(n,1);
E_struct  = zeros(n,1);
E_ratio   = zeros(n,1);
qc_pass   = false(n,1);

T_start    = zeros(n,1);
T_end      = zeros(n,1);
mu_end     = zeros(n,1);
clamp_hits = zeros(n,1);

Dp_mm_col   = repmat(Utils.getfield_default(params,'Dp_mm',NaN), n,1);
mu_ref_col  = repmat(Utils.getfield_default(params,'mu_ref',NaN), n,1);

worstPFA = -inf; worstPFA_name = '';
worstIDR = -inf; worstIDR_name = '';

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
    PFA(k)    = m_nom.PFA_top;
    IDR(k)    = m_nom.IDR_max;
    dP95(k)   = m_nom.dP_orf_q95;
    Qcap95(k) = m_nom.Qcap_ratio_q95;
    cav_pct(k)= m_nom.cav_pct;
    zeta1_hot(k)       = Utils.getfield_default(m_nom,'zeta1_hot',NaN);
    z2_over_z1_hot(k)  = Utils.getfield_default(m_nom,'z2_over_z1_hot',NaN);
    P_mech(k)          = Utils.getfield_default(m_nom,'P_mech',NaN);
    Re_max(k)          = Utils.getfield_default(m_nom,'Re_max',NaN);
    Q_q95(k)  = Utils.getfield_default(m_nom,'Q_q95',NaN);
    Q_q50(k)  = Utils.getfield_default(m_nom,'Q_q50',NaN);
    dP50(k)   = Utils.getfield_default(m_nom,'dP_orf_q50',NaN);
    x10_max_D(k) = Utils.getfield_default(m_nom,'x10_max_D',Utils.getfield_default(m_nom,'x10_pk_D',NaN));
    a10abs_max_D(k) = Utils.getfield_default(m_nom,'a10abs_max_D',Utils.getfield_default(m_nom,'a10abs_pk_D',NaN));
    E_orifice(k) = Utils.getfield_default(m_nom,'E_orifice_full',NaN);
    E_struct(k)  = Utils.getfield_default(m_nom,'E_struct_full',NaN);
    E_ratio(k)   = Utils.getfield_default(m_nom,'E_ratio_full',NaN);
    qc_pass(k)   = out.qc_pass;

    if PFA(k) > worstPFA
        worstPFA = PFA(k);
        worstPFA_name = out.name;
    end
    if IDR(k) > worstIDR
        worstIDR = IDR(k);
        worstIDR_name = out.name;
    end
end

%% Özet Tablo
% Hesaplanan metrikleri tabloya dönüştür
summary = struct();

summary.table = table(names, scale, SaT1, t5, t95, coverage, rank_score, policy_col, order_col, cooldown_col, ...
    PFA, IDR, dP95, Qcap95, cav_pct, zeta1_hot, z2_over_z1_hot, P_mech, Re_max, ...
    Q_q95, Q_q50, dP50, x10_max_D, a10abs_max_D, E_orifice, E_struct, E_ratio, qc_pass, ...
    T_start, T_end, mu_end, clamp_hits, Dp_mm_col, mu_ref_col, ...
    'VariableNames', {'name','scale','SaT1','t5','t95','coverage','rank_score','policy','order','cooldown_s', ...
    'PFA','IDR','dP95','Qcap95','cav_pct','zeta1_hot','z2_over_z1_hot','P_mech','Re_max', ...
    'Q_q95','Q_q50','dP50','x10_max_D','a10abs_max_D','E_orifice','E_struct','E_ratio','qc_pass', ...
    'T_start','T_end','mu_end','clamp_hits','Dp_mm','mu_ref'});

summary.all_out = all_out;
%% QC Kontrolü
% QC eşiklerine göre sonuçların değerlendirilmesi
% --- summary.csv kullanıcıları için QC bayrakları ve sebep kodları ---
thr = opts.thr;
ok_T    = summary.table.T_end   <= thr.T_end_max;
ok_mu   = summary.table.mu_end  >= thr.mu_end_min;
ok_dP   = summary.table.dP95    <= thr.dP95_max;
ok_Qcap = summary.table.Qcap95  <  thr.Qcap95_max;
ok_cav  = summary.table.cav_pct == 0;
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
    fprintf('Worst PFA: %s\n', worstPFA_name);
    fprintf('Worst IDR: %s\n', worstIDR_name);
end

if do_export
    export_results(outdir, scaled, params, opts, summary, all_out);
    if ~Utils.getfield_default(opts,'quiet',false)
        diary off;
    end
end

end
