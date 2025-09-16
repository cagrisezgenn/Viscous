function [paths, env] = init_environment(rootDir)
%INIT_ENVIRONMENT Prepare paths and shared configuration struct.
%   [PATHS, ENV] = INIT_ENVIRONMENT(ROOTDIR) returns a cell array of
%   directories under ROOTDIR that may need to be added to the MATLAB
%   path, as well as a struct ENV containing fields used to share data
%   between functions.
%
%   The function does not modify the global MATLAB path or save it.

if nargin < 1 || isempty(rootDir)
    rootDir = pwd;
end

% Generate path list without altering MATLAB path
p = genpath(rootDir);
paths = regexp(p, pathsep, 'split');
paths = paths(~cellfun('isempty', paths));

% Shut down any existing parallel pool without warning if Parallel
% Toolbox is available.
try
    delete(gcp('nocreate'));
catch
end

% Shared environment struct
env = struct();
env.LOG = struct('verbose_decode', false);
env.PARETO = struct('J1',[],'J2',[],'F',[],'Pen',[], ...
    'set',[],'x',{{}},'feas',[]);
env.PSA_PREP_CACHE = [];
end
