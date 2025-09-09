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
%   OPTS.mu_factors ve OPTS.mu_weights parametreleri RUN_BATCH_WINDOWED
%   fonksiyonundaki varsayılanlarla aynıdır. OPTS.rng_seed "random"
%   sıralamasının tekrarlanabilirliğini sağlar.

%% Girdi ve Varsayılan Ayarlar
if nargin < 3, opts = struct(); end
if ~isfield(opts,'mu_factors'), opts.mu_factors = [0.75 1.00 1.25]; end
if ~isfield(opts,'mu_weights'), opts.mu_weights = [0.2 0.6 0.2]; end
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
    imode = Utils.getfield_default(opts,'IM_mode','');
    band  = Utils.getfield_default(opts,'band_fac',[NaN NaN]);
    sb    = Utils.getfield_default(opts,'s_bounds',[NaN NaN]);
    if ~quiet
        fprintf('Run @ %s | outdir=%s\n', ts, outdir);
        fprintf('IM_mode=%s, band=[%.3g,%.3g], s_bounds=[%.2f,%.2f]\n', imode, band(1), band(2), sb(1), sb(2));
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
    tg = Utils.getfield_default(params,'toggle_gain',NaN);
    tgv = tg(:); tgmin = min(tgv); tgmed = median(tgv); tgmax = max(tgv);
    if ~quiet
        fprintf('Hydraulics: n_orf=%g, d_o~=%g m, A_o=%s, Qcap_big=%g, hA=%g, resFactor=%g, toggle_gain[min/med/max]=[%g %g %g]\n', ...
            n_orf, d_o, mat2str(size(A_o)), Qcap_big, hA, resFactor, tgmin, tgmed, tgmax);
    end
catch ME
    warning('run_policies_windowed echo params: %s', ME.message);
end

nRec = numel(scaled);

%% Temel Koşu
% Delta hesapları ve worst\_first sıralaması için baz koşu
base_opts = opts; base_opts.thermal_reset = 'each'; base_opts.order = 'natural';
base_opts.do_export = false;
[base_summary, base_all] = run_batch_windowed(scaled, params, base_opts);
basePFA_w_mean = mean(base_summary.table.PFA_w);
baseIDR_w_mean = mean(base_summary.table.IDR_w);
baseTend_worst = max(base_summary.table.T_end_worst);
base_dP95_worst_max = max(base_summary.table.dP95_worst);
base_Qcap95_worst_max = max(base_summary.table.Qcap95_worst);
base_cav_worst_max = max(base_summary.table.cav_pct_worst);
base_mu_end_worst_min = min(base_summary.table.mu_end_worst);
base_qc_pass = sum(base_summary.table.qc_all_mu);
base_qc_n = height(base_summary.table);
if ~quiet
    fprintf(['BASE (each/natural): PFA_w=%.3g, IDR_w=%.3g, dP95_worst=%.3g MPa, Qcap95_worst=%.2f, ' ...
            'cav%%_worst=%.1f, T_end_worst=%.1f C, mu_end_worst=%.2f, qc_rate=%d/%d\n'], ...
        basePFA_w_mean, baseIDR_w_mean, base_dP95_worst_max/1e6, base_Qcap95_worst_max, base_cav_worst_max*100, baseTend_worst,
        base_mu_end_worst_min, base_qc_pass, base_qc_n);
end

%% Sıra Ön Hesapları
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
        case 'pfa_w'
            rk = cellfun(@(s) s.weighted.PFA_top, base_all);
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

%% Kombinasyon Döngüsü
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

            qc.pass_fraction = mean(summary.table.qc_all_mu);
            qc.n = height(summary.table);

            % Kombinasyon için özetler
            curPFA_w_mean = mean(summary.table.PFA_w);
            curIDR_w_mean = mean(summary.table.IDR_w);
            curTend_worst = max(summary.table.T_end_worst);
            deltas = struct('PFA_w', curPFA_w_mean - basePFA_w_mean, ...
                            'IDR_w', curIDR_w_mean - baseIDR_w_mean, ...
                            'T_end_worst', curTend_worst - baseTend_worst);

            % En kötü durumların kaydı
            [worstPFA, idxP] = max(summary.table.PFA_worst);
            nP = summary.table.name{idxP};
            muP = summary.table.which_mu_PFA(idxP);
            if ~quiet, fprintf('Worst PFA (%s,%s,mu=%.2f): %s\n', pol, ord, muP, nP); end
            [worstIDR, idxI] = max(summary.table.IDR_worst); %#ok<NASGU>
            nI = summary.table.name{idxI};
            muI = summary.table.which_mu_IDR(idxI);
            if ~quiet, fprintf('Worst IDR (%s,%s,mu=%.2f): %s\n', pol, ord, muI, nI); end

            % Politika karşılaştırma logu
            tolPFA = 0.15 * basePFA_w_mean;
            tolIDR = 0.15 * baseIDR_w_mean;
            passPFA = abs(deltas.PFA_w) <= tolPFA;
            passIDR = abs(deltas.IDR_w) <= tolIDR;
            if ~quiet
                fprintf('Delta vs base: dPFA_w=%.4g (|d|<=%.4g? %d), dIDR_w=%.4g (|d|<=%.4g? %d), qc_rate=%.2f\n', ...
                    deltas.PFA_w, tolPFA, passPFA, deltas.IDR_w, tolIDR, passIDR, qc.pass_fraction);
            end

            % Tek satırlık özet
            n_pass = sum(summary.table.qc_all_mu);
            n_tot  = height(summary.table);
            PFAw   = curPFA_w_mean; IDRw = curIDR_w_mean; T_end_worst_max = curTend_worst;
            pass_flag_15pct = 'OK'; if ~(passPFA && passIDR), pass_flag_15pct = 'FAIL'; end
            if ~quiet
                fprintf(['policy=%s | order=%s | cd=%ds | PFA_w=%.3g (d=%+.2f%%) | IDR_w=%.3g (d=%+.2f%%) | ' ...
                        'T_end_worst=%.1f C | qc_rate=%d/%d %s\n'], ...
                    pol, ord, cdval, PFAw, 100*(PFAw-basePFA_w_mean)/max(basePFA_w_mean,eps), ...
                    IDRw, 100*(IDRw-baseIDR_w_mean)/max(baseIDR_w_mean,eps), T_end_worst_max, n_pass, n_tot, pass_flag_15pct);
                % En kötü iki kaydı bildir
                try
                    kshow = min(2, height(summary.table));
                    [~,ix] = maxk(summary.table.PFA_worst, kshow);
                    if kshow==2
                        fprintf('  worst2(PFA): %s | %s\n', summary.table.name{ix(1)}, summary.table.name{ix(2)});
                    elseif kshow==1
                        fprintf('  worst2(PFA): %s | -\n', summary.table.name{ix(1)});
                    end
                catch ME
                    warning('run_policies_windowed worst2: %s', ME.message);
                end
            end

            % Opsiyonel clamp özeti
            try
                if ~quiet
                    tot_clamps = sum(summary.table.clamp_hits);
                    nz_idx = find(summary.table.clamp_hits>0);
                    nz_names = summary.table.name(nz_idx);
                    if ~isempty(nz_names)
                        fprintf('clamp_hits total=%d | records: %s\n', tot_clamps, strjoin(nz_names.', ', '));
                    else
                        fprintf('clamp_hits total=%d\n', tot_clamps);
                    end
                end
            catch ME
                warning('run_policies_windowed clamp summary: %s', ME.message);
            end

            P(end+1) = struct('policy',pol,'order',ord,'cooldown_s',cdval, ...
                'summary',summary.table,'qc',qc,'deltas',deltas); %#ok<AGROW>
        end
    end
end

%% Sonuçların Kaydedilmesi
if do_export
    export_results(outdir, scaled, params, opts, base_summary, base_all, P);
end
if ~quiet, fprintf('Saved to %s\n', outdir); end
if do_export && ~quiet, diary off; end
end
