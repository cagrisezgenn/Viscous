function metr = compute_metrics_windowed(t, x, a_rel, ag, diag, story_height, win)
%COMPUTE_METRICS_WINDOWED Pencereli metrikleri hesaplar.
%   Metrikler hem pencere içinde hem de tüm kayıtta hesaplanır.
%
%   Girişler:
%     t, x, a_rel, ag  - zaman ve çözümler
%     diag             - mck_with_damper çıktısı
%     story_height     - kat yüksekliği (m)
%     win              - make_arias_window çıktısı
%
%   Çıkış:
%     metr yapısı: PFA_top, IDR_max, çeşitli akış/enerji/termal metrikler.

idx = win.idx;
if numel(t) ~= size(x,1)
    error('Zaman ve çözüm boyutları uyumsuz.');
end

% PFA üst kat
if isempty(a_rel)
    PFA_top = NaN;
else
    a_top = a_rel(:,end) + ag(:);
    PFA_top = max(abs(a_top(idx)));
end

% IDR maksimumu
if isempty(x)
    IDR_max = NaN;
else
    drift = diff(x,1,2);
    IDR_max = max(max(abs(drift(idx,:))))/max(story_height,eps);
end

metr = struct('PFA_top',PFA_top,'IDR_max',IDR_max);

% --- diag tabanlı metrikler (varsa) ---
q95 = @(v) prctile(v,95);

fields = {'dP_orf','Q','Qcap_ratio','cav','story_force'};
for k = 1:numel(fields)
    f = fields{k};
    if isfield(diag,f)
        val = diag.(f);
        if isvector(val), val = val(:); end
        metr.([f '_q95']) = q95(abs(val(idx,:)));
        if strcmp(f,'cav')
            metr.cav_pct = 100*mean(val(idx,:)<0);
        end
    else
        metr.([f '_q95']) = NaN;
        if strcmp(f,'cav')
            metr.cav_pct = NaN;
        end
    end
end

% Enerji
if isfield(diag,'E_orifice')
    E_orf_series = diag.E_orifice;
    E_orifice_full = E_orf_series(end);
    E_orifice_win  = E_orf_series(find(idx,1,'last')) - E_orf_series(find(idx,1,'first'));
else
    E_orifice_win = NaN; E_orifice_full = NaN;
end

if isfield(diag,'E_struct')
    E_struct_series = diag.E_struct;
    E_struct_full = E_struct_series(end);
    E_struct_win  = E_struct_series(find(idx,1,'last')) - E_struct_series(find(idx,1,'first'));
else
    E_struct_win = NaN; E_struct_full = NaN;
end

E_win = E_orifice_win + E_struct_win;
E_full = E_orifice_full + E_struct_full;
E_ratio = E_orifice_full / max(E_struct_full,eps);
E_cov = E_win / max(E_full,eps);

metr.E_orifice = E_orifice_win;
metr.E_struct  = E_struct_win;
metr.E_ratio   = E_ratio;
metr.E_win_over_full = E_cov;

% Termal
if isfield(diag,'T_oil')
    metr.T_oil_end = diag.T_oil(find(idx,1,'last'));
else
    metr.T_oil_end = NaN;
end
if isfield(diag,'mu')
    metr.mu_end = diag.mu(find(idx,1,'last'));
else
    metr.mu_end = NaN;
end

% Modal metrikler
if isfield(diag,'zeta1_hot')
    metr.zeta1_hot = diag.zeta1_hot;
else
    metr.zeta1_hot = NaN;
end
if isfield(diag,'zeta2_over_zeta1_hot')
    metr.z2_over_z1 = diag.zeta2_over_zeta1_hot;
else
    metr.z2_over_z1 = NaN;
end

% Kapsama (E_win/E_full) zaten hesaplandı
end
