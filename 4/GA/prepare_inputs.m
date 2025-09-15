function [scaled, params, T1] = prepare_inputs(optsGA)
%PREPARE_INPUTS Çalışma alanı bağımsız giriş hazırlığı.
%   [SCALED, PARAMS, T1] = PREPARE_INPUTS(OPTSGA) gerekli yolları ekler,
%   parametre yapısını oluşturur ve yer hareketi verilerini yükler.
%   OPTSGA.load_opts varsa load_ground_motions'a iletilir.

if nargin < 1 || isempty(optsGA)
    optsGA = struct;
end

% Gerekli yollar için opsiyonel setup çağrısı
if exist('setup','file') == 2
    setup;
end

% Parametreler ve T1 hesapla (script çağrısı fonksiyon kapsamı içinde)
parametreler;  %#ok<NOPRT>

% Yer hareketi verisini ölçekle
load_opts = Utils.getfield_default(optsGA, 'load_opts', []);
if isempty(load_opts)
    [~, scaled] = load_ground_motions(T1);
else
    [~, scaled] = load_ground_motions(T1, load_opts);
end

end
