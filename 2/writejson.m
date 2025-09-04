function writejson(data, filename)
%WRITEJSON Minimal JSON writer using JSONENCODE.
try
    txt = jsonencode(data);
    fid = fopen(filename,'w');
    if fid~=-1
        fwrite(fid, txt);
        fclose(fid);
    end
catch
    % silently ignore
end
end
