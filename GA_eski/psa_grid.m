function [tPSA, agPSA] = psa_grid(t, ag, down_dt, warn_tail)
%PSA_GRID Optional resampling for PSA (downsample only; no upsample).
%   [tPSA, agPSA] = PSA_GRID(t, ag, down_dt) resamples the acceleration
%   record AG defined at times T onto a coarser grid with spacing DOWN_DT.
%   If DOWN_DT is empty or not greater than the native time step, the input
%   vectors are returned unchanged. Resampling uses PCHIP interpolation.
%   If WARN_TAIL is true, emit a warning when DOWN_DT does not evenly divide
%   the record length, resulting in truncation.

    t = t(:); ag = ag(:);
    if nargin < 3 || isempty(down_dt) || down_dt <= 0
        tPSA = t; agPSA = ag; return;
    end
    if nargin < 4 || isempty(warn_tail)
        warn_tail = false;
    end
    dt0 = median(diff(t), 'omitnan');
    if down_dt <= dt0 * (1 + 1e-12)
        tPSA = t; agPSA = ag; return;
    end
    tPSA  = (t(1):down_dt:t(end)).';
    % Clip to input domain to avoid extrapolation
    tPSA(tPSA < t(1)) = t(1);
    tPSA(tPSA > t(end)) = t(end);
    agPSA = interp1(t, ag, tPSA, 'pchip');
    if warn_tail
        tail = t(end) - tPSA(end);
        if tail > max(eps(t(end)), eps(down_dt))
            warning('psa_grid:tailTruncated', ...
                'down_dt does not evenly divide record length; last %.6g s truncated', tail);
        end
    end
end
