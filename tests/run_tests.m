function run_tests
% Basit doğrulama testleri

% 1. Arias penceresi testi
Fs = 100; t = (0:1/Fs:10).'; ag = sin(2*pi*1*t);
win = make_arias_window(t, ag);
assert(win.coverage > 0.9, 'coverage düşük');

% 2. Metrik hesaplama testi
x = [0.01*sin(2*pi*1*t), 0.02*sin(2*pi*1*t)];
a_rel = x;
diag = struct('E_orifice',cumsum(0.1*ones(size(t))), ...
              'E_struct',cumsum(0.2*ones(size(t))), ...
              'T_oil',20+0.1*t, 'mu',0.9*ones(size(t)));
metr = compute_metrics_windowed(t, x, a_rel, ag, diag, 3.0, win);
assert(isfield(metr,'PFA_top'));

fprintf('run_tests: tüm testler geçti.\n');
end
