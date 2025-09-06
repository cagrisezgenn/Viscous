classdef Utils
%UTILS Collect small helper functions as static methods.
    methods(Static)
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
    end
end

