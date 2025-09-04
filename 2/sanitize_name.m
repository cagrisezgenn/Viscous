function s = sanitize_name(s)
%SANITIZE_NAME Replace non-alphanumeric characters for safe file names.
%   S2 = SANITIZE_NAME(S) replaces characters outside [a-zA-Z0-9_\- ] with '_'.
if ~ischar(s) && ~isstring(s)
    s = char(s);
end
s = regexprep(char(s),'[^a-zA-Z0-9_\- ]','_');
end
