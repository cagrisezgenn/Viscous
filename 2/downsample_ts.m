function ts_ds = downsample_ts(ts, ds)
%DOWNSAMPLE_TS Downsample all numeric fields in a struct.
%   TS_DS = DOWNSAMPLE_TS(TS, DS) keeps every DS-th row of numeric array fields.
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
