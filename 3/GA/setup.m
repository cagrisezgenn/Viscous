function oldPath = setup()
% SETUP proje klasörlerini MATLAB yoluna ekler.
% Test ve çalıştırma için gerekli dizinleri dahil eder.

    % Project root is the folder containing this setup.m
    projectRoot = fileparts(mfilename('fullpath'));

    %% Yolların Oluşturulması
    % Build a genpath and filter out unwanted folders
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
        % Exclude common non-code directories and examples
        if strcmpi(name, 'out') || strcmpi(name, '.git') || startsWith(name, '.') || strcmpi(name, 'examples')
            keep(i) = false; %#ok<AGROW>
        end
    end

    filteredPath = strjoin(parts(keep), pathsep);

    %% Yolun Eklenmesi
    % Return old path and add the filtered one
    oldPath = path;
    addpath(filteredPath);
end
