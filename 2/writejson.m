function writejson(fname, data)
%WRITEJSON Write MATLAB data to a JSON file.
%   WRITEJSON(FNAME, DATA) encodes DATA as JSON and writes it to file FNAME.

text = jsonencode(data);
fid = fopen(fname,'w');
if fid == -1, error('writejson:open','Cannot open %s', fname); end
fwrite(fid, text);
fclose(fid);

