function [xbest, fbest, output, pop, scores, exitflag] = ga_call_compat(fhandle, lb, ub, IntCon, opts)
%GA_CALL_COMPAT Wrapper to call GA if available; otherwise fallback search.
%  Tries multiple GA signatures for toolbox/version compatibility.
%  If GA is unavailable, performs a lightweight random/LHS search using the
%  provided bounds and options (PopulationSize/MaxGenerations/InitialPop).

    % Defaults
    if nargin < 5 || isempty(opts), opts = struct(); end
    if nargin < 4, IntCon = []; end

    % Try the official GA if available
    if exist('ga','file') == 2
        try
            % Newer signatures (6 outputs)
            [xbest, fbest, exitflag, output, pop, scores] = ...
                ga(fhandle, numel(lb), [], [], [], [], lb, ub, [], IntCon, opts);
            return;
        catch
            try
                % Mid signatures (4 outputs)
                [xbest, fbest, exitflag, output] = ...
                    ga(fhandle, numel(lb), [], [], [], [], lb, ub, [], IntCon, opts);
                pop = []; scores = [];
                return;
            catch
                try
                    % Old minimal signatures (2 outputs)
                    [xbest, fbest] = ga(fhandle, numel(lb), [], [], [], [], lb, ub, [], IntCon, opts);
                    exitflag = []; output = struct(); pop = []; scores = [];
                    return;
                catch
                    % Fall through to manual fallback
                end
            end
        end
    end

    % Manual fallback: simple LHS/random search
    Npop = get_opt(opts, 'PopulationSize', 16);
    G    = get_opt(opts, 'MaxGenerations', 8);
    P0   = get_opt(opts, 'InitialPopulationMatrix', []);

    % Build initial population
    if ~isempty(P0)
        P = ensure_pop(P0, lb, ub, Npop);
    else
        P = lhs_population(lb, ub, Npop);
    end

    xbest = []; fbest = inf;
    pop = []; scores = [];

    for g = 1:G
        fvals = nan(size(P,1),1);
        for i = 1:size(P,1)
            fvals(i) = fhandle(P(i,:));
            if fvals(i) < fbest
                fbest = fvals(i); xbest = P(i,:);
            end
        end
        pop = P; scores = fvals;
        % Next generation: around current best + random LHS
        P = jitter_around_best(xbest, lb, ub, Npop, max(0.05, 0.5*(G-g)/G));
    end

    output = struct('message','ga_call_compat: fallback_random_search');
    exitflag = NaN;
end

function val = get_opt(opts, name, def)
    try
        if isstruct(opts)
            if isfield(opts,name), val = opts.(name); else, val = def; end
        else
            % optim.options.GaOptions or similar
            if isprop(opts,name) || ismethod(opts,name)
                val = opts.(name);
            else
                val = def;
            end
        end
    catch
        val = def;
    end
end

function P = lhs_population(lb,ub,N)
    D = numel(lb); P = zeros(N,D);
    for d=1:D
        edges = linspace(0,1,N+1);
        centers = (edges(1:end-1)+edges(2:end))/2;
        P(:,d) = lb(d) + centers(randperm(N))' .* (ub(d)-lb(d));
    end
end

function P = ensure_pop(P0, lb, ub, N)
    % Clamp and pad/trim to N rows
    P = max(lb(:)'.*ones(size(P0,1),1), min(ub(:)'.*ones(size(P0,1),1), P0));
    if size(P,1) < N
        P = [P; lhs_population(lb,ub,N-size(P,1))]; %#ok<AGROW>
    elseif size(P,1) > N
        P = P(1:N,:);
    end
end

function P = jitter_around_best(xb, lb, ub, N, frac)
    if isempty(xb)
        P = lhs_population(lb,ub,N); return;
    end
    D = numel(xb); P = zeros(N,D);
    span = (ub - lb) .* frac;
    for i=1:N
        xi = xb + (rand(1,D)-0.5).*2.*span;
        P(i,:) = max(lb(:)', min(ub(:)', xi));
    end
end


