function h = md5(str)
%MD5 Compute MD5 hash of input string, returning lowercase hex digest.
%   h = MD5(str) returns a 32-char lowercase hexadecimal MD5 digest.

    % Normalize input to character row vector
    if isstring(str)
        str = char(str);
    elseif isnumeric(str)
        str = char(str);
    end
    if ischar(str)
        str = str(:).';
    else
        error('md5:InvalidInput', 'Input must be char or string.');
    end

    % Prefer Octave's hash if available; otherwise, use Java in MATLAB
    if exist('OCTAVE_VERSION','builtin') ~= 0
        try
            % Octave: default output is hex already
            h = lower(hash('md5', str));
            return;
        catch
            % Fall through to Java implementation
        end
    end

    % Java MessageDigest (works in MATLAB; available with JVM)
    md = java.security.MessageDigest.getInstance('MD5');
    % One-shot digest call with MATLAB uint8 -> Java byte[] conversion
    if isempty(str)
        bytes = md.digest();
    else
        % Use Java to create UTF-8 bytes from the MATLAB char string
        jstr = java.lang.String(str);
        jbytes = jstr.getBytes('UTF-8');
        bytes = md.digest(jbytes);
    end
    % Convert Java byte[] to hex using Java BigInteger to avoid MATLAB conversion issues
    bi = java.math.BigInteger(1, bytes);
    h = char(bi.toString(16));
    if numel(h) < 32
        h = [repmat('0', 1, 32 - numel(h)) h];
    end
    h = lower(h);
end

