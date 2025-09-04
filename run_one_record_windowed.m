function out = run_one_record_windowed(rec, params, opts)
%RUN_ONE_RECORD_WINDOWED Tek deprem kaydını çalıştırır (pencereli).
%   rec    : scaled(k) elemanı (t, ag, name, scale, IM)
%   params : parametreler.m çıktıları (M,K,C0, ...)
%   opts   : yapı (window, sim, thermal_reset, ...)
%
%   out alanları: name, scale, SaT1, win, metr, metr0, diag, flags,
%   which_peak.

if nargin<3, opts = struct(); end

% --- Varsayılan seçenekler ---
default.window.mode = 'scaled';
default.window.pad  = 0.5;
default.sim.mode    = 'full';
default.thermal_reset = 'each';
default.cooldown_s    = 0;
opts = merge_struct(default, opts);

% 1) Arias penceresi
win = make_arias_window(rec.t, rec.ag, 'pad', opts.window.pad, 'mode', opts.window.mode);

% 2) Dampersiz çözüm
[x0, a_rel0] = lin_MCK(rec.t, rec.ag, params.M, params.C0, params.K);
metr0.PFA_top = max(abs(a_rel0(win.idx,end) + rec.ag(win.idx)));
drift0 = diff(x0,1,2);
metr0.IDR_max = max(max(abs(drift0(win.idx,:))))/params.story_height;

% 3) Damperli çözüm
use_orf = getfield(params,'use_orifice',true);
use_thermal = getfield(params,'use_thermal',true);
[x, a_rel, diag] = mck_with_damper(rec.t, rec.ag, params.M, params.C0, params.K, ...
    params.k_sd, params.c_lam0, use_orf, params.orf, params.rho, params.Ap, params.A_o, params.Qcap_big, params.mu_ref, ...
    use_thermal, params.thermal, params.T0_C, params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, params.Lgap, ...
    params.cp_oil, params.cp_steel, params.steel_to_oil_mass_ratio, params.toggle_gain, params.story_mask, ...
    params.n_dampers_per_story, params.resFactor, params.cfg);

% 4) Metrikler
metr = compute_metrics_windowed(rec.t, x, a_rel, rec.ag, diag, params.story_height, win);

% 5) Hangi zamanda tepe?
t_win = rec.t(win.idx);
[~,k_pfa] = max(abs(a_rel(win.idx,end) + rec.ag(win.idx)));
which_peak.PFA.time = t_win(k_pfa);
which_peak.PFA.value = metr.PFA_top;

% IDR tepe
if ~isempty(x)
    drift = diff(x(win.idx,:),1,2);
    [mx,idxmx] = max(abs(drift(:)));
    [row,col] = ind2sub(size(drift),idxmx);
    which_peak.IDR.time  = t_win(row);
    which_peak.IDR.story = col;
    which_peak.IDR.value = mx/params.story_height;
else
    which_peak.IDR = struct('time',NaN,'story',NaN,'value',NaN);
end

% 6) Tutarlılık kontrolü (a_abs ≈ a_rel + ag)
if ~isempty(x)
    r = randi(numel(t_win));
    idx_glob = find(win.idx);
    k = idx_glob(r);
    dt = mean(diff(rec.t));
    v_top = gradient(x(:,end), dt);
    a_from_x = gradient(v_top, dt);
    err = abs((a_from_x(k) + rec.ag(k)) - (a_rel(k,end) + rec.ag(k)));
    if err > 1e-6
        warning('a_abs kontrolü başarısız: %.3e', err);
    end
end

out = struct();
out.name   = rec.name;
out.scale  = rec.scale;
out.SaT1   = rec.IM;
out.win    = win;
out.metr   = metr;
out.metr0  = metr0;
out.diag   = diag;
out.flags  = struct('low_arias',win.flag_low_arias);
out.which_peak = which_peak;
end

function s = merge_struct(a,b)
    s = a;
    f = fieldnames(b);
    for i=1:numel(f), s.(f{i}) = b.(f{i}); end
end
