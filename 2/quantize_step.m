function y = quantize_step(x, step)
%QUANTIZE_STEP Snap scalar/vector x to a uniform grid with spacing STEP.
%   y = step * round(x/step)
    if nargin < 2 || isempty(step)
        y = x;
        return;
    end
    y = step * round(x ./ step);
end

