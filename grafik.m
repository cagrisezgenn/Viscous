%% ================================================================
%  Sonuç Grafikleri ve Kısa Özet
%  Bu betik, <code>damperlinon.m</code> tarafından çağrılır. Çözümlerin
%  karşılaştırma grafiğini üretir ve maksimum yer değiştirmeleri yazdırır.
%  Aşağıdaki değişkenlerin ana betikte tanımlanmış olması beklenir:
%    t, x10_0, x10_lin, x10_orf, T1, t5, t95, use_thermal
%% ================================================================

%% Grafik: 10. kat yer değiştirme eğrileri
figure('Name','10. Kat yer değiştirme — ham ivme (ODE-only)','Color','w');
plot(t, x10_0 ,'k','LineWidth',1.4); hold on;
plot(t, x10_lin,'b','LineWidth',1.1);
plot(t, x10_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('x_{10}(t) [m]');
title(sprintf('10-Kat | T1=%.3f s | Arias [%.3f, %.3f] s', T1, t5, t95));
legend('Dampersiz','Lineer damper','Orifisli damper','Location','best');

%% Kısa özet metni
fprintf('x10_max  (dampersiz)   = %.4g m\n', max(abs(x10_0)));
fprintf('x10_max  (lineer)      = %.4g m\n', max(abs(x10_lin)));
if use_thermal
    suffix = '+termal';
else
    suffix = '';
end
fprintf('x10_max  (orifisli%s)  = %.4g m\n', suffix, max(abs(x10_orf)));if exist('diag','var') && isstruct(diag) && isfield(diag,'dT_est')
    fprintf('Termal döngü: ΔT_est=%.2f K | c_lam(final)=%.3e N·s/m\n', diag.dT_est, diag.c_lam);
end
fprintf('Not: orifis modelini kapatmak için use_orifice=false; termali kapatmak için use_thermal=false.\n');

