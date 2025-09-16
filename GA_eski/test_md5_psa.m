% Minimal runtime checks for md5 and psa_hash
addpath(pwd);

% Known MD5 vectors
a = md5("");
b = md5("abc");
fprintf("md5('')   = %s\n", a);
fprintf("md5('abc')= %s\n", b);
assert(strcmp(a, "d41d8cd98f00b204e9800998ecf8427e"));
assert(strcmp(b, "900150983cd24fb0d6963f7d28e17f72"));

% psa_hash quick check
k = psa_hash([1 2], [3 4], [0.1 0.2], 0.05);
fprintf("psa_hash sample = %s\n", k);
assert(ischar(k) && numel(k) == 32);
assert(all(ismember(k, '0123456789abcdef')));

disp('All tests passed.');


