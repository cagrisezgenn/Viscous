function s = sanitize_name(name)
%SANITIZE_NAME Replace non-alphanumeric characters for safe filenames.
%   S = SANITIZE_NAME(NAME) replaces characters outside a-z, A-Z, 0-9, underscore,
%   dash and space with underscores for use in file and directory names.

s = regexprep(name, '[^a-zA-Z0-9_\- ]', '_');

