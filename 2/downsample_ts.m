function ts_ds = downsample_ts(ts, ds)
%DOWNSAMPLE_TS Downsample all numeric fields of a time-series struct.
%   TS_DS = DOWNSAMPLE_TS(TS, DS) selects every DS-th sample along the first
%   dimension for all numeric vector or matrix fields in TS.  DS defaults to 5.

if nargin < 2 || isempty(ds), ds = 5; end
fields = fieldnames(ts);
ts_ds = struct();
for i = 1:numel(fields)
    f = fields{i};
    val = ts.(f);
    if isnumeric(val) && ~isscalar(val)
        sz = size(val);
        if sz(1) == 1
            ts_ds.(f) = val(:,1:ds:end);
        else
            ts_ds.(f) = val(1:ds:end,:);
        end
    else
        ts_ds.(f) = val;
    end
end

