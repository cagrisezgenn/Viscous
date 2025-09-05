function oldPath = setup()
%SETUP Add project folders to MATLAB path for testing/running.
%   oldPath = SETUP() adds the project root and subfolders (excluding
%   common output/hidden folders) to the MATLAB path and returns the
%   previous path so you can restore it with path(oldPath) when finished.

    % Project root is the folder containing this setup.m
    projectRoot = fileparts(mfilename('fullpath'));

    % Build a genpath and filter out unwanted folders
    rawPath = genpath(projectRoot);
    parts = strsplit(rawPath, pathsep);

    keep = true(size(parts));
    for i = 1:numel(parts)
        folder = parts{i};
        if isempty(folder)
            keep(i) = false; %#ok<AGROW>
            continue
        end
        [~, name] = fileparts(folder);
        % Exclude common non-code directories
        if strcmpi(name, 'out') || strcmpi(name, '.git') || startsWith(name, '.')
            keep(i) = false; %#ok<AGROW>
        end
    end

    filteredPath = strjoin(parts(keep), pathsep);

    % Return old path and add the filtered one
    oldPath = path;
    addpath(filteredPath);
end

