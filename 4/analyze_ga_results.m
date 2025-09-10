function R = analyze_ga_results(outdir)
% ANALYZE_GA_RESULTS Son GA klasörünü okuyup hızlı karar desteği üretir.
%   R = ANALYZE_GA_RESULTS()  -> out/ga_* içinden en yeni klasörü alır
%   R = ANALYZE_GA_RESULTS(outdir)
% Çıktı R yapılandırma önerileri ve temel istatistikleri içerir.

if nargin < 1 || isempty(outdir)
    dd = dir(fullfile('out','ga_*'));
    assert(~isempty(dd), 'No GA outputs under out/ga_*');
    [~,ix] = max([dd.datenum]);
    outdir = fullfile(dd(ix).folder, dd(ix).name);
end

frontCsv = fullfile(outdir,'ga_front.csv');
assert(exist(frontCsv,'file')==2, 'ga_front.csv not found in %s', outdir);
T = readtable(frontCsv);

% Baseline ilk satır olabilir; analizden hariç tut
if height(T) >= 2 && any(isnan(T.f1(1)))
    T = T(2:end,:);
end

% Eşikler: snapshot varsa oradan al; yoksa Utils default
thr = Utils.default_qc_thresholds(struct());
snap = fullfile(outdir,'..','snapshot.mat'); % parent may not match; try local too
snap1 = fullfile(outdir,'snapshot.mat');
try
    if exist(snap1,'file')==2
        S = load(snap1);
    elseif exist(snap,'file')==2
        S = load(snap);
    else
        S = struct();
    end
    if isfield(S,'thr'), thr = Utils.default_qc_thresholds(S.thr); end
catch
end

% İhlal bayrakları
bad_dP   = T.dP95_worst   > thr.dP95_max;
bad_Qcap = T.Qcap95_worst >= thr.Qcap95_max;
bad_cav  = T.cav_pct_worst > max(thr.cav_pct_max, 0);
bad_T    = T.T_end_worst  > thr.T_end_max;
bad_mu   = T.mu_end_worst < thr.mu_end_min;

n = height(T);
rate = @(x) mean(x(:));
rates = struct('dP',rate(bad_dP), 'Qcap',rate(bad_Qcap), 'cav',rate(bad_cav), ...
               'T',rate(bad_T), 'mu',rate(bad_mu));

% Hedef fonksiyonlar
f1m = mean(T.f1); f2m = mean(T.f2);

% Basit öneri mantığı
Wrec = struct('dP',1,'Qcap',1,'cav',3,'T',1,'mu',0.5);
recs = {};
if rates.cav > 0.15
    recs{end+1} = 'Increase penalty_weights.cav (e.g., 3→4 or 5).';
end
if rates.Qcap > 0.15
    recs{end+1} = 'Increase penalty_weights.Qcap and consider raising hA_W_perK or lowering PF_gain.';
end
if rates.dP > 0.15
    recs{end+1} = 'Increase penalty_weights.dP and consider narrowing CdInf upper bound.';
end
if rates.T > 0.15
    recs{end+1} = 'Increase penalty_weights.T and consider higher hA_W_perK bound.';
end
if rates.mu > 0.15
    recs{end+1} = 'Increase penalty_weights.mu or shift mu_factors to [0.7 1.0 1.4].';
end
if isempty(recs)
    recs{end+1} = 'Keep current weights; consider local refine on knee/topK.';
end

R = struct();
R.outdir = outdir;
R.f1_mean = f1m; R.f2_mean = f2m;
R.violation_rates = rates;
R.penalty_suggestion = Wrec;
R.recommendations = recs;

fprintf('Analysis of %s\n', outdir);
fprintf('  mean f1=%.4g, mean f2=%.4g\n', f1m, f2m);
fprintf('  violation rates: dP=%.0f%%, Qcap=%.0f%%, cav=%.0f%%, T=%.0f%%, mu=%.0f%%\n', ...
    100*rates.dP, 100*rates.Qcap, 100*rates.cav, 100*rates.T, 100*rates.mu);
for i=1:numel(recs)
    fprintf('  - %s\n', recs{i});
end

end

