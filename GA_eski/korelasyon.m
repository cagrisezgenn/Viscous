% KORELASYON.M
% trial_log.csv dosyasini okuyup, ilk 30 degisken (set1..set5)
% ile secilen metrikler arasindaki Pearson korelasyonunu hesaplar.
% Esik uygulanmaz; gecerli (sonlu) tum korelasyonlar korelasyon.csv'ye yazilir.

% Dosya konumu: oncelikle out/trial_log.csv, yoksa trial_log.csv
candidates = {fullfile('out','trial_log.csv'), 'trial_log.csv'};
fp = '';
for k = 1:numel(candidates)
    if exist(candidates{k}, 'file')
        fp = candidates{k};
        break;
    end
end
if isempty(fp)
    error('trial_log.csv bulunamadi (out/trial_log.csv veya trial_log.csv).');
end

fprintf('Okunuyor: %s\n', fp);
T = readtable(fp);

% Ilk 30 sutun set degiskenleri, 31..end metrikler
nvars = width(T);
if nvars < 30
    error('trial_log.csv en az 30 set sutunu icermeli.');
end

idx_set = 1:30;
varnames = T.Properties.VariableNames;

% Sadece bu metrikler: acc3_max_abs, x10_max_abs
targets = {'acc3_max_abs','x10_max_abs'};
idx_met = [];
met_names = {};
for ii = 1:numel(targets)
    k = find(strcmp(varnames, targets{ii}), 1, 'first');
    if ~isempty(k)
        idx_met(end+1) = k; %#ok<AGROW>
        met_names{end+1} = varnames{k}; %#ok<AGROW>
    end
end
if isempty(idx_met)
    warning('Hedef metrikler bulunamadi: %s', strjoin(targets, ', '));
end

set_names = varnames(idx_set);

rows = {};
for i = 1:numel(idx_set)
    xi = T{:, idx_set(i)};
    for j = 1:numel(idx_met)
        yj = T{:, idx_met(j)};
        mask = isfinite(xi) & isfinite(yj);
        n = sum(mask);
        if n >= 3 && std(xi(mask)) > 0 && std(yj(mask)) > 0
            C = corrcoef(xi(mask), yj(mask));
            r = C(1,2);
        else
            r = NaN;
        end
        if isfinite(r)
            rows(end+1, :) = {set_names{i}, met_names{j}, r, n, abs(r)}; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    fprintf('Gecerli (sonlu) korelasyon bulunamadi.\n');
    % Bos tablo da yazilabilir
    Out = cell2table(cell(0,4), 'VariableNames', {'set_var','metric','corr','N'});
else
    % |corr| azalan sirada
    Rtab = cell2table(rows, 'VariableNames', {'set_var','metric','corr','N','abs_corr'});
    Rtab = sortrows(Rtab, 'abs_corr', 'descend');
    Rtab.abs_corr = [];
    Out = Rtab;
end

out_fp = 'korelasyon.csv';
writetable(Out, out_fp);
fprintf('Yazildi: %s (satir sayisi=%d)\n', out_fp, height(Out));

