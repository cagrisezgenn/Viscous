function oldPath = setup()
% SETUP proje klasörlerini MATLAB yoluna ekler.
% Test ve çalıştırma için gerekli dizinleri dahil eder.

    % Proje kökü, bu setup.m dosyasını içeren klasördür
    projectRoot = fileparts(mfilename('fullpath'));

    %% Yolların Oluşturulması
    % Bir genpath oluştur ve istenmeyen klasörleri filtrele
    rawPath = genpath(projectRoot);
    parts = strsplit(rawPath, pathsep);

    %% İstenmeyen Klasörlerin Filtrelenmesi
    keep = true(size(parts));
    for i = 1:numel(parts)
        folder = parts{i};
        if isempty(folder)
            keep(i) = false; %#ok<AGROW>
            continue
        end
        [~, name] = fileparts(folder);
        % Kod içermeyen yaygın dizinleri ve örnekleri hariç tut
        if strcmpi(name, 'out') || strcmpi(name, '.git') || startsWith(name, '.') || strcmpi(name, 'examples')
            keep(i) = false; %#ok<AGROW>
        end
    end

    filteredPath = strjoin(parts(keep), pathsep);

    %% Yolun Eklenmesi
    % Eski yolu döndür ve filtrelenmiş olanı ekle
    oldPath = path;
    addpath(filteredPath);
end
