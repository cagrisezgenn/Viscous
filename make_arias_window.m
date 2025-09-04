function win = make_arias_window(t, ag, varargin)
%MAKE_ARIAS_WINDOW Arias yoğunluk penceresi hesaplar.
%   win = MAKE_ARIAS_WINDOW(t, ag) varsayılan olarak Arias %5-%95
%   penceresini ve 0.5 s tamponu kullanır. Ek parametreler
%   'p1','p2','pad','mode' anahtarları ile geçilebilir.
%
%   Çıktı alanları:
%     t5, t95      - Arias yoğunluğunun %p1 ve %p2 zamanları
%     pad          - kullanılan tampon (s)
%     t_start      - pencere başlangıcı
%     t_end        - pencere sonu
%     idx          - pencereyi gösteren mantıksal vektör
%     coverage     - a_g^2 integralinin kapsanan oranı
%     flag_low_arias - toplam enerji çok düşükse true
%
%   Kenar durumları için:
%     * Eğer kayıt kısa ise (t95-t5 < 5 s) pad=0.25 s yapılır.
%     * Arias toplamı ~0 ise flag_low_arias=true döner.

p = inputParser;
p.addParameter('p1',0.05,@(x)isscalar(x)&&x>=0&&x<=1);
p.addParameter('p2',0.95,@(x)isscalar(x)&&x>=0&&x<=1);
p.addParameter('pad',0.5,@(x)isscalar(x)&&x>=0);
p.addParameter('mode','scaled',@(x)ischar(x)||isstring(x));
p.parse(varargin{:});
pr = p.Results;

[t5, t95] = arias_win(t, ag, pr.p1, pr.p2);
pad = pr.pad;
if t95 - t5 < 5
    pad = 0.25; % kısa kayıt
end

% t_start ve t_end hesapla
if isempty(t)
    t_start = 0; t_end = 0; idx = false(size(t));
else
    t_start = max(t(1), t5 - pad);
    t_end   = min(t(end), t95 + pad);
    idx = (t >= t_start) & (t <= t_end);
end

% kapsama oranı
E_total = trapz(t, ag.^2);
if E_total <= eps
    coverage = 0;
    flag_low_arias = true;
else
    coverage = trapz(t(idx), ag(idx).^2) / E_total;
    flag_low_arias = (E_total < 1e-4);
end

win = struct('t5',t5,'t95',t95,'pad',pad,...
    't_start',t_start,'t_end',t_end,'idx',idx,...
    'coverage',coverage,'flag_low_arias',flag_low_arias);
end
