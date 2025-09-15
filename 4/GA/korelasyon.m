function S = korelasyon(outdir, varargin)
% KORELASYON  GA çıktılarında değişken ↔ çıktı ilişkilerini hızlı analiz eder.
%   S = KORELASYON(outdir, ...) 4/GA/out altındaki ga_*/ga_front.csv dosyalarını
%   toplayıp temel korelasyonlar, kısmi korelasyonlar, lineer/ansambıl önem
%   analizleri ve (mevcutsa) kısmi bağımlılık/SHAP benzeri açıklamalar üretir.
%
%   Girdiler (opsiyonel):
%     outdir: Varsayılan '4/GA/out'.
%     İsim-değer çiftleri:
%       'Responses'     : Hücre dizi; hedef çıktı kolon adları
%       'Predictors'    : Hücre dizi; tahminci kolon adları
%       'SaveDir'       : Sonuç figürleri için klasör (vars: fullfile(outdir,'corr'))
%       'UseEconometrics': true ise corrplot kullanılır (varsa)
%       'Spearman'      : true ise Spearman korelasyon da hesaplanır
%       'Verbose'       : true ise özet çıktı konsola yazılır
%
%   Çıktı (S):
%     .T        : Birleştirilmiş tablo
%     .R_pearson: Pearson korelasyon matrisi (predictors+responses)
%     .R_spearman: Spearman korelasyon matrisi (opsiyonel)
%     .partial  : Yapısal (predictor | diğer predictor’lar) kısmi kor. (yanıt başına)
%     .lm       : Standartlaştırılmış lineer modeller (yanıt başına)
%     .ens      : Ansambıl modeller ve predictor önemleri (yanıt başına)
%     .figs     : Kaydedilen figür dosya yolları
%
% Notlar:
% - Statistics and Machine Learning Toolbox fonksiyonları (corr, partialcorr,
%   fitlm, fitrensemble, predictorImportance) tercih edilir.
% - Econometrics Toolbox ‘corrplot’ varsa kullanılabilir; yoksa ısı haritası çizilir.

% Varsayılanlar ve seçenekler
if nargin < 1 || isempty(outdir)
    % Betiğin bulunduğu klasörü referans al: <repo>/4/GA → <repo>/4/GA/out
    here = fileparts(mfilename('fullpath'));
    outdir = fullfile(here, 'out');
end
p = inputParser;
p.addParameter('Responses', { ...
    'f1','f2','PFA_mean','IDR_mean','T_end','mu_end', ...
    'dP95','Qcap95','PF_p95','cav_pct', ...
    'x10_max_damperli','a10abs_max_damperli', ...
    'Q_q50','Q_q95','dP50', ...
    'energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio','P_mech_sum' ...
});
p.addParameter('Predictors', { ...
    'd_o_mm','n_orf','PF_tau','PF_gain','PF_t_on', ...
    'Cd0','CdInf','p_exp','Lori_mm','hA_W_perK', ...
    'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref' ...
});
p.addParameter('SaveDir', fullfile(outdir,'corr'));
p.addParameter('UseEconometrics', true, @(x)islogical(x)||isnumeric(x));
p.addParameter('Spearman', true, @(x)islogical(x)||isnumeric(x));
p.addParameter('Verbose', true, @(x)islogical(x)||isnumeric(x));
p.parse(varargin{:});
opt = p.Results; opt.UseEconometrics = logical(opt.UseEconometrics);
opt.Spearman = logical(opt.Spearman); opt.Verbose = logical(opt.Verbose);

if ~exist(opt.SaveDir,'dir'), mkdir(opt.SaveDir); end

% 1) ga_front.csv dosyalarını topla
% Öncelik: outdir/ga_front.csv (tekleştirilmiş derleme) → alt klasörler/ga_*/ga_front.csv
T_all = table();
singleFile = fullfile(outdir, 'ga_front.csv');
if exist(singleFile,'file') == 2
    try
        T = readtable(singleFile);
        T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);
        T_all = T;
    catch ME
        warning('Tek dosya okunamadı (%s): %s', singleFile, ME.message);
    end
end
if isempty(T_all)
    files = dir(fullfile(outdir, 'ga_*', 'ga_front.csv'));
    assert(~isempty(files), 'korelasyon: ga_front.csv bulunamadı: %s', outdir);
    for k = 1:numel(files)
        try
            T = readtable(fullfile(files(k).folder, files(k).name));
            % Kolon adlarını temizle (boşluk/simge toleransı için)
            T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);
            T_all = [T_all; T]; %#ok<AGROW>
        catch ME
            warning('Dosya okunamadı (%s): %s', files(k).name, ME.message);
        end
    end
end
assert(~isempty(T_all), 'korelasyon: Okunabilir tablo bulunamadı.');

% 2) Kullanılabilir kolonları kesişimle seç
pred_all = intersect(opt.Predictors, T_all.Properties.VariableNames, 'stable');
resp_all = intersect(opt.Responses,  T_all.Properties.VariableNames, 'stable');
assert(~isempty(pred_all) && ~isempty(resp_all), 'korelasyon: predictor/response kolonları bulunamadı.');

% Sayısal tablo oluştur (NaN toleransı)
Xraw = T_all(:, pred_all);  Yraw = T_all(:, resp_all);
for c = 1:width(Xraw)
    if ~isnumeric(Xraw{:,c}), Xraw{:,c} = double(string(Xraw{:,c})); end
end
for c = 1:width(Yraw)
    if ~isnumeric(Yraw{:,c}), Yraw{:,c} = double(string(Yraw{:,c})); end
end
X = Xraw{:,:}; Y = Yraw{:,:};
% NaN satırları at
keep = all(isfinite(X),2) & all(isfinite(Y),2);
X = X(keep,:); Y = Y(keep,:);
T_use = [array2table(X,'VariableNames',pred_all), array2table(Y,'VariableNames',resp_all)];

% 3) Korelasyonlar
[R_pearson, P_pearson] = corrcoef(T_use{:,:}, 'Rows','pairwise');
R_spear = []; P_spear = [];
if opt.Spearman
    try
        [R_spear,P_spear] = corr(T_use{:,:}, 'Type','Spearman','Rows','pairwise');
    catch
        warning('Spearman hesaplanamadı (corr yok?)');
    end
end

% 4) Kısmi korelasyon (yanıt başına; diğer predictor’lara koşullu)
partial = struct();
try
    for j = 1:numel(resp_all)
        y = Y(:,j);
        pc = nan(1, numel(pred_all)); pv = pc;
        for i = 1:numel(pred_all)
            Xi = X(:,i); Xc = X; Xc(:,i) = [];
            [pc(i), pv(i)] = partialcorr(Xi, y, Xc, 'Rows','pairwise');
        end
        partial.(resp_all{j}) = table(pred_all', pc', pv', 'VariableNames', {'predictor','partial_r','pvalue'});
    end
catch ME
    warning('partialcorr başarısız: %s', ME.message);
end

% 5) Standartlaştırılmış lineer modeller (beta ~ etki büyüklüğü)
lm = struct();
try
    Xz = (X - mean(X,1)) ./ max(std(X,[],1), eps);
    for j = 1:numel(resp_all)
        yz = (Y(:,j) - mean(Y(:,j))) / max(std(Y(:,j)), eps);
        mdl = fitlm(Xz, yz, 'Intercept', true, 'VarNames', [pred_all, resp_all(j)]);
        lm.(resp_all{j}) = mdl;
    end
catch ME
    warning('fitlm başarısız: %s', ME.message);
end

% 6) Ansambıl modeller (ağaç/Bagging) ve önem
ens = struct();
try
    for j = 1:numel(resp_all)
        yj = Y(:,j);
        mdlE = fitrensemble(X, yj, 'Method','Bag', 'NumLearningCycles', 200, ...
            'Learners', templateTree('MinLeafSize', max(3,round(size(X,1)/200))), ...
            'KFold', 5);
        % Out-of-fold modelden önem
        try
            imp = kfoldPredictorImportance(mdlE);
        catch
            % Eski sürümlerde kfoldPredictorImportance yoksa tek model ile deneyelim
            mdls = fitrensemble(X, yj, 'Method','Bag', 'NumLearningCycles', 200);
            imp = predictorImportance(mdls);
        end
        ens.(resp_all{j}) = struct('model', mdlE, 'importance', imp(:), 'predictors', {pred_all});
    end
catch ME
    warning('fitrensemble başarısız: %s', ME.message);
end

% 7) Görselleştirme ve kayıt
figs = struct();
try
    vnames = T_use.Properties.VariableNames;
    if opt.UseEconometrics && exist('corrplot','file') == 2
        % corrplot için yorumlayıcıları kapat (none) ve ham isimleri kullan
        oldTI = get(groot,'defaultTextInterpreter');
        oldLI = get(groot,'defaultLegendInterpreter');
        % ColorbarInterpreter tüm sürümlerde mevcut olmayabilir
        try, oldCI = get(groot,'defaultColorbarInterpreter'); hasCI = true; catch, hasCI = false; oldCI = ''; end
        set(groot,'defaultTextInterpreter','none');
        set(groot,'defaultLegendInterpreter','none');
        if hasCI, try, set(groot,'defaultColorbarInterpreter','none'); end, end
        f = figure('Name','corrplot');
        corrplot(T_use, 'VarNames', vnames);
        figs.corrplot = fullfile(opt.SaveDir,'corrplot.png');
        saveas(f, figs.corrplot); close(f);
        % eski ayarlara dön
        set(groot,'defaultTextInterpreter',oldTI);
        set(groot,'defaultLegendInterpreter',oldLI);
        if hasCI, try, set(groot,'defaultColorbarInterpreter',oldCI); end, end
    else
        f = figure('Name','corr_heatmap');
        imagesc(R_pearson); colorbar; axis square;
        xticks(1:numel(vnames)); yticks(1:numel(vnames));
        xticklabels(vnames); yticklabels(vnames); xtickangle(45);
        ax = gca; ax.TickLabelInterpreter = 'none';
        title('Pearson Correlation');
        figs.corr_heatmap = fullfile(opt.SaveDir,'corr_heatmap.png');
        saveas(f, figs.corr_heatmap); close(f);
    end
catch ME
    warning('Korelasyon görselleştirme başarısız: %s', ME.message);
end

% Önem grafikleri (yanıt başına ilk 10)
try
    fn = fieldnames(ens);
    for k = 1:numel(fn)
        e = ens.(fn{k});
        if ~isfield(e,'importance') || isempty(e.importance), continue; end
        [~,ord] = sort(e.importance,'descend');
        top = min(10, numel(ord));
        f = figure('Name',['importance_' fn{k}]);
        bar(e.importance(ord(1:top))); grid on;
        xticklabels(e.predictors(ord(1:top))); xtickangle(30);
        ax = gca; ax.TickLabelInterpreter = 'none';
        ylabel('Importance'); title(sprintf('Predictor Importance -> %s', fn{k}));
        figs.(['imp_' fn{k}]) = fullfile(opt.SaveDir, sprintf('imp_%s.png', fn{k}));
        saveas(f, figs.(['imp_' fn{k}])); close(f);
    end
catch ME
    warning('Önem grafikleri başarısız: %s', ME.message);
end

% 8) Özet yazdır
if opt.Verbose
    fprintf('[korelasyon] Kayıt sayısı = %d | Pred=%d | Resp=%d\n', size(T_use,1), numel(pred_all), numel(resp_all));
    try
        fprintf('Pearson korelasyon (ilk 5x5):\n'); disp(R_pearson(1:min(end,5),1:min(end,5)));
    catch
    end
    try
        yname = resp_all{1};
        if isfield(lm,yname)
            fprintf('fitlm(%s) özeti:\n', yname);
            disp(lm.(yname).Coefficients);
        end
    catch
    end
end

% 9) Çıktı paketleme
S = struct('T', T_use, 'R_pearson', R_pearson, 'P_pearson', P_pearson, ...
    'R_spearman', R_spear, 'P_spearman', P_spear, ...
    'partial', partial, 'lm', lm, 'ens', ens, 'figs', figs);

%% 10) CSV çıktılarını yaz
try
    % Veri birleşimi
    data_csv = fullfile(opt.SaveDir, 'data_merged.csv');
    writetable(T_use, data_csv);

    % Korelasyon/p-değer tabloları (Pearson)
    names_all = T_use.Properties.VariableNames;
    TcorrP = local_mat_to_table(R_pearson, names_all);
    TpvalP = local_mat_to_table(P_pearson, names_all);
    writetable(TcorrP, fullfile(opt.SaveDir, 'corr_pearson.csv'));
    writetable(TpvalP, fullfile(opt.SaveDir, 'pval_pearson.csv'));

    % Spearman varsa
    if ~isempty(R_spear)
        TcorrS = local_mat_to_table(R_spear, names_all);
        TpvalS = local_mat_to_table(P_spear, names_all);
        writetable(TcorrS, fullfile(opt.SaveDir, 'corr_spearman.csv'));
        writetable(TpvalS, fullfile(opt.SaveDir, 'pval_spearman.csv'));
    end

    % Kısmi korelasyon (yanıt başına)
    fns = fieldnames(partial);
    for i = 1:numel(fns)
        Tpar = partial.(fns{i});
        writetable(Tpar, fullfile(opt.SaveDir, sprintf('partial_%s.csv', fns{i})));
    end

    % Lineer modeller (standartlaştırılmış)
    fns = fieldnames(lm);
    for i = 1:numel(fns)
        mdl = lm.(fns{i});
        C = mdl.Coefficients;
        % RowNames → sütun olarak ekle
        try
            terms = C.Properties.RowNames;
        catch
            terms = strcat('term', string((1:height(C))'));
        end
        Tcoef = C; %#ok<NASGU>
        Tcoef = addvars(C, terms, 'Before',1, 'NewVariableNames','term');
        writetable(Tcoef, fullfile(opt.SaveDir, sprintf('lm_%s.csv', fns{i})));
    end

    % Ansambıl önem (yanıt başına)
    fns = fieldnames(ens);
    for i = 1:numel(fns)
        e = ens.(fns{i});
        if isfield(e,'importance') && ~isempty(e.importance)
            Timp = table(e.predictors(:), e.importance(:), 'VariableNames', {'predictor','importance'});
            writetable(Timp, fullfile(opt.SaveDir, sprintf('ens_%s.csv', fns{i})));
        end
    end
catch ME
    warning('CSV yazma hatası: %s', ME.message);
end

end

%% === Yerel yardımcılar ===
function Tm = local_mat_to_table(M, names)
    % Kare matris + isim listesi → bir CSV dostu tablo
    Tm = array2table(M, 'VariableNames', names);
    Tm = addvars(Tm, names(:), 'Before', 1, 'NewVariableNames', 'name');
end
