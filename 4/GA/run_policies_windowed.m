function P = run_policies_windowed(scaled, params, opts)
%RUN_POLICIES_WINDOWED Termal reset politikalarını ve kayıt sıralarını değerlendirir.
%   P = RUN_POLICIES_WINDOWED(SCALED, PARAMS, OPTS) termal reset
%   politikaları ile kayıt işleme sıralarının kombinasyonlarını çalıştırır.
%   OPTS.policies {'each','carry','cooldown'} alt kümesini, OPTS.orders ise
%   {'natural','random','worst_first'} seçeneklerini içerir. 'cooldown'
%   politikasında OPTS.cooldown_s_list bekleme sürelerini belirtir. Her
%   kombinasyon için politika, sıra, bekleme süresi, özet tablo ve QC
%   istatistikleri döndürülür.
%
%   OPTS.rng_seed "random" sıralamasının tekrarlanabilirliğini sağlar.

%% Girdi ve Varsayılan Ayarlar
if nargin < 3, opts = struct(); end
if ~isfield(opts,'policies'), opts.policies = {'each','carry','cooldown'}; end
if ~isfield(opts,'orders'), opts.orders = {'natural','random','worst_first'}; end
if ~isfield(opts,'cooldown_s_list'), opts.cooldown_s_list = [60 180 300]; end
if ~isfield(opts,'rng_seed'), opts.rng_seed = 42; end
if ~isfield(opts,'rank_metric'), opts.rank_metric = 'E_orifice_win'; end

% QC eşikleri (Utils ile varsayılanlara tamamlanır)
if ~isfield(opts,'thr'), opts.thr = struct(); end
opts.thr = Utils.default_qc_thresholds(opts.thr);

% Sessizlik/çıktı bayrakları ve çıktı dizini
quiet = isfield(opts,'quiet') && opts.quiet;
% Politika koşuları için varsayılan olarak sonuçlar dışa aktarılır
do_export = ~isfield(opts,'do_export') || opts.do_export;
ts = datestr(now,'yyyymmdd_HHMMSS_FFF'); outdir = fullfile('out', ts);
if do_export && ~quiet
    if ~exist(outdir,'dir'), mkdir(outdir); end
    diary(fullfile(outdir,'console.log'));
else
    if ~exist(outdir,'dir'), mkdir(outdir); end
end

% Konsol başlığı
try
    if ~quiet
        fprintf('Run @ %s | outdir=%s\n', ts, outdir);
        fprintf('policies=%s | orders=%s | cooldown_s_list=%s | rng_seed=%d\n', ...
            strjoin(opts.policies,','), strjoin(opts.orders,','), sprintf('%d ', opts.cooldown_s_list), opts.rng_seed);
        if isfield(opts,'TRIM_names') && ~isempty(opts.TRIM_names)
            fprintf('TRIM: %s\n', strjoin(opts.TRIM_names,', '));
        end
        fprintf('QC thr: dP95<=%.1f MPa, Qcap95<%.2f, cav%%=%g, T_end<=%g C, mu_end>=%0.2f Pa*s\n', ...
            opts.thr.dP95_max/1e6, opts.thr.Qcap95_max, opts.thr.cav_pct_max*100, opts.thr.T_end_max, opts.thr.mu_end_min);
    end
catch ME
    warning('run_policies_windowed header: %s', ME.message);
end

% Hidrolik ve termal temel parametreleri yazdır
try
    n_orf = NaN; if isfield(params,'A_o'), n_orf = numel(params.A_o); end
    A_o = Utils.getfield_default(params,'A_o',NaN);
    d_o = NaN; try, d_o = sqrt(4*mean(A_o)/pi); catch, end
    Qcap_big = Utils.getfield_default(params,'Qcap_big',NaN);
    hA = NaN; if isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), hA = params.thermal.hA_W_perK; end
    resFactor = Utils.getfield_default(params,'resFactor',NaN);
    if ~quiet
        fprintf('Hydraulics: n_orf=%g, d_o~=%g m, A_o=%s, Qcap_big=%g, hA=%g, resFactor=%g\n', n_orf, d_o, mat2str(size(A_o)), Qcap_big, hA, resFactor);
    end
catch ME
    warning('run_policies_windowed echo params: %s', ME.message);
end

%% Temel Koşu
% Delta hesapları ve worst_first sıralaması için baz koşu
base_opts = opts; base_opts.thermal_reset = 'each'; base_opts.order = 'natural';
base_opts.do_export = false;
[base_summary, base_all] = run_batch_windowed(scaled, params, base_opts);
basePFA_mean = mean(base_summary.table.PFA);
baseIDR_mean = mean(base_summary.table.IDR);
baseTend_max = max(base_summary.table.T_end);
base_dP95_max = max(base_summary.table.dP95);
base_Qcap95_max = max(base_summary.table.Qcap95);
base_cav_max = max(base_summary.table.cav_pct);
base_mu_end_min = min(base_summary.table.mu_end);
base_qc_pass = sum(base_summary.table.qc_pass);
base_qc_n = height(base_summary.table);
if ~quiet
    fprintf(['BASE (each/natural): PFA=%.3g, IDR=%.3g, dP95=%.3g MPa, Qcap95=%.2f, ' ...
            'cav%%=%.1f, T_end_max=%.1f C, mu_end_min=%.2f, qc_rate=%d/%d\n'], ...
        basePFA_mean, baseIDR_mean, base_dP95_max/1e6, base_Qcap95_max, base_cav_max*100, baseTend_max, ...
        base_mu_end_min, base_qc_pass, base_qc_n);
end

orders_struct = compute_orders(opts, scaled, base_all, quiet);
base_metrics = struct('PFA_mean', basePFA_mean, 'IDR_mean', baseIDR_mean, 'T_end_max', baseTend_max);
P = run_combinations(scaled, params, opts, orders_struct, base_metrics, quiet);

%% Sonuçların Kaydedilmesi
if do_export
    export_results(outdir, scaled, params, opts, base_summary, base_all, P);
end
if ~quiet, fprintf('Saved to %s\n', outdir); end
if do_export && ~quiet, diary off; end
end

function orders_struct = compute_orders(opts, scaled, base_all, quiet)
%COMPUTE_ORDERS Farklı sıralama stratejilerini hazırlar
nRec = numel(scaled);
orders_struct.natural = 1:nRec;
if any(strcmp(opts.orders,'random'))
    rng(opts.rng_seed);
    orders_struct.random = randperm(nRec);
end
if any(strcmp(opts.orders,'worst_first'))
    switch lower(opts.rank_metric)
        case 'pfa_top'
            rk = cellfun(@(s) s.metr.PFA_top, base_all);
        case 'idr_max'
            rk = cellfun(@(s) s.metr.IDR_max, base_all);
        otherwise
            rk = cellfun(@(s) s.metr.E_orifice_win, base_all);
    end
    [~,idx] = sort(rk,'descend');
    orders_struct.worst_first = idx;
    try
        rank_names = {scaled(idx).name};
        if ~quiet
            fprintf('worst_first ranking by %s: %s\n', opts.rank_metric, strjoin(rank_names, ', '));
        end
    catch ME
        warning('run_policies_windowed ranking: %s', ME.message);
    end
end
end

function P = run_combinations(scaled, params, opts, orders_struct, base_metrics, quiet)
%RUN_COMBINATIONS Politika-sıra kombinasyonlarını değerlendirir
P = struct('policy',{},'order',{},'cooldown_s',{},'summary',{},'qc',{},'deltas',{});
for ip = 1:numel(opts.policies)
    pol = opts.policies{ip};
    for io = 1:numel(opts.orders)
        ord = opts.orders{io};
        if strcmp(pol,'cooldown')
            cds = opts.cooldown_s_list(:)';
        else
            cds = NaN;
        end
        for ic = 1:numel(cds)
            cdval = cds(ic);
            perm = orders_struct.(ord);
            scaled_run = scaled(perm);
            run_opts = opts;
            if isfield(run_opts,'cooldown_s'), run_opts = rmfield(run_opts,'cooldown_s'); end
            run_opts.order = ord;
            run_opts.thermal_reset = pol;
            run_opts.do_export = false;
            if strcmp(pol,'cooldown'), run_opts.cooldown_s = cdval; end
            [summary, ~] = run_batch_windowed(scaled_run, params, run_opts);

            qc.pass_fraction = mean(summary.table.qc_pass);
            qc.n = height(summary.table);

            curPFA_mean = mean(summary.table.PFA);
            curIDR_mean = mean(summary.table.IDR);
            curTend_max = max(summary.table.T_end);
            deltas = struct('PFA', curPFA_mean - base_metrics.PFA_mean, ...
                            'IDR', curIDR_mean - base_metrics.IDR_mean, ...
                            'T_end', curTend_max - base_metrics.T_end_max);

            report_combination(pol, ord, cdval, summary, deltas, qc, base_metrics, quiet);

            P(end+1) = struct('policy',pol,'order',ord,'cooldown_s',cdval, ...
                'summary',summary.table,'qc',qc,'deltas',deltas); %#ok<AGROW>
        end
    end
end
end

function report_combination(pol, ord, cdval, summary, deltas, qc, base_metrics, quiet)
%REPORT_COMBINATION Kombinasyon sonuçlarını loglar
if quiet, return; end
[~, idxP] = max(summary.table.PFA);
fprintf('Worst PFA (%s,%s): %s\n', pol, ord, summary.table.name{idxP});
[~, idxI] = max(summary.table.IDR);
fprintf('Worst IDR (%s,%s): %s\n', pol, ord, summary.table.name{idxI});

% Politika karşılaştırma logu
tolPFA = 0.15 * base_metrics.PFA_mean;
tolIDR = 0.15 * base_metrics.IDR_mean;
passPFA = abs(deltas.PFA) <= tolPFA;
passIDR = abs(deltas.IDR) <= tolIDR;
fprintf('Delta vs base: dPFA=%.4g (|d|<=%.4g? %d), dIDR=%.4g (|d|<=%.4g? %d), qc_rate=%.2f\n', ...
    deltas.PFA, tolPFA, passPFA, deltas.IDR, tolIDR, passIDR, qc.pass_fraction);

% Tek satırlık özet
n_pass = sum(summary.table.qc_pass);
n_tot  = height(summary.table);
PFAm   = mean(summary.table.PFA);
IDRm   = mean(summary.table.IDR);
T_end_max = max(summary.table.T_end);
pass_flag_15pct = 'OK'; if ~(passPFA && passIDR), pass_flag_15pct = 'FAIL'; end
fprintf(['policy=%s | order=%s | cd=%ds | PFA=%.3g (d=%+.2f%%) | IDR=%.3g (d=%+.2f%%) | ' ...
        'T_end_max=%.1f C | qc_rate=%d/%d %s\n'], ...
    pol, ord, cdval, PFAm, 100*(PFAm-base_metrics.PFA_mean)/max(base_metrics.PFA_mean,eps), ...
    IDRm, 100*(IDRm-base_metrics.IDR_mean)/max(base_metrics.IDR_mean,eps), T_end_max, n_pass, n_tot, pass_flag_15pct);
try
    kshow = min(2, height(summary.table));
    [~,ix] = maxk(summary.table.PFA, kshow);
    if kshow==2
        fprintf('  worst2(PFA): %s | %s\n', summary.table.name{ix(1)}, summary.table.name{ix(2)});
    elseif kshow==1
        fprintf('  worst2(PFA): %s | -\n', summary.table.name{ix(1)});
    end
catch ME
    warning('run_policies_windowed worst2: %s', ME.message);
end
try
    tot_clamps = sum(summary.table.clamp_hits);
    nz_idx = find(summary.table.clamp_hits>0);
    nz_names = summary.table.name(nz_idx);
    if ~isempty(nz_names)
        fprintf('clamp_hits total=%d | records: %s\n', tot_clamps, strjoin(nz_names.', ', '));
    else
        fprintf('clamp_hits total=%d\n', tot_clamps);
    end
catch ME
    warning('run_policies_windowed clamp summary: %s', ME.message);
end
end
