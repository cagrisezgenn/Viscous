classdef Utils
%UTILS Collect small helper functions as static methods.
    methods(Static)
        %% Recompute derived damper parameters (Phase 1 helper)
        function params = recompute_damper_params(params)
            % Recompute Ap, k_p, k_sd, c_lam0 based on basic inputs.
            % Accepts optional mm-suffixed fields and converts to meters.
            if ~isstruct(params), return; end

            % mm -> m conversions if present
            if isfield(params,'Dp_mm'),    params.Dp    = params.Dp_mm/1000; end
            if isfield(params,'d_w_mm'),   params.d_w   = params.d_w_mm/1000; end
            if isfield(params,'D_m_mm'),   params.D_m   = params.D_m_mm/1000; end
            if isfield(params,'Lori_mm'),  params.Lori  = params.Lori_mm/1000; end
            if isfield(params,'orf') && isfield(params.orf,'d_o_mm')
                params.orf.d_o = params.orf.d_o_mm/1000;
            end

            req = {'Dp','d_w','D_m','n_turn','mu_ref','Lori','Lgap','Kd','Ebody','Gsh'};
            if ~all(isfield(params,req)) || ~isfield(params,'orf') || ~isfield(params.orf,'d_o')
                return; % missing fields; skip
            end

            Ap = pi * params.Dp^2 / 4;
            k_h = params.Kd * Ap^2 / params.Lgap;
            k_s = params.Ebody * Ap / params.Lgap;
            k_hyd = 1 / (1/max(k_h,eps) + 1/max(k_s,eps));
            k_p = params.Gsh * params.d_w^4 / (8 * params.n_turn * params.D_m^3);
            k_sd = k_hyd + k_p;
            c_lam0 = 12 * params.mu_ref * params.Lori * Ap^2 / (max(params.orf.d_o,eps)^4);

            params.Ap = Ap;
            params.k_p = k_p;
            params.k_sd = k_sd;
            params.c_lam0 = c_lam0;
        end
        function y = softmin(a,b,epsm)
            y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
        end

        function w = pf_weight(t, cfg)
            w = cfg.on.pressure_force * (1 - exp(-max(t - cfg.PF.t_on, 0) ./ max(cfg.PF.tau, 1e-6)));
        end

        function [x,a] = lin_MCK(t,ag,M,C,K)
            n = size(M,1); r = ones(n,1);
            agf = griddedInterpolant(t,ag,'linear','nearest');
            odef = @(tt,z)[ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - M*r*agf(tt) ) ];
            z0 = zeros(2*n,1);
            opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
            sol = ode15s(odef,[t(1) t(end)],z0,opts);
            z = deval(sol,t).';
            x = z(:,1:n);
            a = ( -(M\(C*z(:,n+1:end).' + K*z(:,1:n).')).' - ag.*r.' );
        end

        function win = make_arias_window(t, ag, varargin)
            p = inputParser;
            p.addParameter('p1',0.05,@(x)isscalar(x) && x>=0 && x<=1);
            p.addParameter('p2',0.95,@(x)isscalar(x) && x>=0 && x<=1);
            p.addParameter('pad',0.5,@(x)isscalar(x) && x>=0);
            p.parse(varargin{:});
            p1 = p.Results.p1; p2 = p.Results.p2; pad = p.Results.pad;
            IA = cumtrapz(t, ag.^2);
            IA_tot = IA(end);
            % Guard against degenerate/very low-energy records
            if ~(isfinite(IA_tot)) || IA_tot <= eps
                t_start = t(1); t_end = t(end);
                idx = true(size(t));
                win = struct('t5',t_start,'t95',t_end,'pad',0, ...
                             't_start',t_start,'t_end',t_end,'idx',idx, ...
                             'coverage',1.0,'flag_low_arias',true);
                return;
            end
            IA_norm = IA / IA_tot;
            [t_unique, iu] = unique(IA_norm);
            interp_t = t(iu);
            t5  = interp1(t_unique, interp_t, p1, 'linear');
            t95 = interp1(t_unique, interp_t, p2, 'linear');
            dur = t95 - t5;
            if dur < 5, pad = 0.25; end
            t_start = max(t(1),  t5  - pad);
            t_end   = min(t(end), t95 + pad);
            idx = (t >= t_start) & (t <= t_end);
            coverage = trapz(t(idx), ag(idx).^2) / IA_tot;
            flag_low_arias = coverage < 0.90;
            win = struct('t5',t5,'t95',t95,'pad',pad, ...
                         't_start',t_start,'t_end',t_end,'idx',idx, ...
                         'coverage',coverage,'flag_low_arias',flag_low_arias);
        end

        function y = quantize_step(x, step)
            if nargin < 2 || isempty(step)
                y = x; return;
            end
            y = step * round(x ./ step);
        end

        function v = getfield_default(S, fname, defaultVal)
            if ~isstruct(S) || ~isfield(S, fname) || isempty(S.(fname))
                v = defaultVal; return;
            end
            val = S.(fname);
            if isnumeric(val)
                if isempty(val) || any(~isfinite(val(:)))
                    v = defaultVal;
                else
                    v = val;
                end
            else
                v = val;
            end
        end

        function writejson(data, filename)
            try
                txt = jsonencode(data);
                fid = fopen(filename,'w');
                if fid~=-1
                    fwrite(fid, txt);
                    fclose(fid);
                end
            catch
            end
        end

        function s2 = sanitize_name(s)
            if ~ischar(s) && ~isstring(s)
                s = char(s);
            end
            s2 = regexprep(char(s),'[^a-zA-Z0-9_\- ]','_');
        end

        function ts_ds = downsample_ts(ts, ds)
            if nargin < 2 || isempty(ds), ds = 5; end
            fns = fieldnames(ts);
            ts_ds = struct();
            for i = 1:numel(fns)
                val = ts.(fns{i});
                if isnumeric(val) && ~isscalar(val)
                    idx = 1:ds:size(val,1);
                    ts_ds.(fns{i}) = val(idx,:,:);
                else
                    ts_ds.(fns{i}) = val;
                end
            end
        end

        function s = tern(c,a,b)
            if c, s=a; else, s=b; end
        end

        function P = initial_pop_grid(lb, ub, N, steps)
            %INITIAL_POP_GRID Create an initial population aligned to variable grids.
            %   P = Utils.initial_pop_grid(lb, ub, N, steps) returns an NxD matrix
            %   where each column i is sampled on a grid with spacing steps(i)
            %   (if steps(i) is NaN, sample uniformly in [lb(i), ub(i)]).
            d = numel(lb);
            P = zeros(N, d);
            for i = 1:d
                if ~isnan(steps(i))
                    if steps(i) <= 0
                        gs = max(ub(i)-lb(i), eps)/10;
                    else
                        gs = steps(i);
                    end
                    grid = lb(i):gs:ub(i);
                    if numel(grid) < 2
                        k = max(2, ceil((ub(i)-lb(i))/max(gs, eps)));
                        grid = linspace(lb(i), ub(i), k);
                    end
                    idx = randi(numel(grid), [N,1]);
                    vals = grid(idx);
                else
                    vals = lb(i) + rand(N,1).*(ub(i)-lb(i));
                end
                P(:,i) = vals;
            end
        end

        %% Default QC Thresholds
        function thr = default_qc_thresholds(optsThr)
            % Return default quality-control thresholds and fill missing fields.
            % Example: thr = Utils.default_qc_thresholds(struct('dP95_max',40e6));

            if nargin < 1 || isempty(optsThr)
                optsThr = struct();
            end

            thr = struct();
            thr.dP95_max   = Utils.getfield_default(optsThr,'dP95_max',50e6);
            thr.Qcap95_max = Utils.getfield_default(optsThr,'Qcap95_max',0.5);
            thr.cav_pct_max= Utils.getfield_default(optsThr,'cav_pct_max',0);
            thr.T_end_max  = Utils.getfield_default(optsThr,'T_end_max',75);
            thr.mu_end_min = Utils.getfield_default(optsThr,'mu_end_min',0.5);

            % Preserve any additional fields
            extra = setdiff(fieldnames(optsThr), fieldnames(thr));
            for ii = 1:numel(extra)
                thr.(extra{ii}) = optsThr.(extra{ii});
            end
        end
    end
end
