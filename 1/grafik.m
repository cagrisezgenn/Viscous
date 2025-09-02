%% ================================================================
%  Sonuç Grafikleri ve Kısa Özet
%  Bu betik, <code>damperlinon.m</code> tarafından çağrılır. Çözümlerin
%  karşılaştırma grafiğini üretir ve maksimum yer değiştirmeleri yazdırır.
%  Aşağıdaki değişkenlerin ana betikte tanımlanmış olması beklenir:
%    t, x10_0, x10_lin, x10_orf, T1, t5, t95, use_thermal

%% Grafik: 10. kat yer değiştirme eğrileri
figure('Name','10. Kat yer değiştirme — ham ivme (ODE-only)','Color','w');
plot(t, x10_0 ,'k','LineWidth',1.4); hold on;
plot(t, x10_lin,'b','LineWidth',1.1);
plot(t, x10_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('x10(t) [m]');
title(sprintf('10-Kat | T1=%.3f s | Arias [%.3f, %.3f] s', T1, t5, t95));
legend('Dampersiz','Lineer damper','Orifisli damper','Location','best');

%% Grafik: 10. kat mutlak ivme
figure('Name','10. Kat mutlak ivme','Color','w');
plot(t, a10_0 ,'k','LineWidth',1.4); hold on;
plot(t, a10_lin,'b','LineWidth',1.1);
plot(t, a10_orf,'r','LineWidth',1.0);
yl = ylim; plot([t5 t5],yl,'k--','HandleVisibility','off');
plot([t95 t95],yl,'k--','HandleVisibility','off');
grid on; xlabel('t [s]'); ylabel('a10abs(t) [m/s^2]');
legend('Dampersiz','Lineer damper','Orifisli damper','Location','best');

%% Grafik: Maksimum göreli kat ötelemeleri (IDR)
drift0    = x0(:,2:end)   - x0(:,1:end-1);
drift_lin = x_lin(:,2:end) - x_lin(:,1:end-1);
drift_orf = x_orf(:,2:end) - x_orf(:,1:end-1);
IDR0      = max(abs(drift0))./story_height;
IDR_lin   = max(abs(drift_lin))./story_height;
IDR_orf   = max(abs(drift_orf))./story_height;
story_ids = 1:(n-1);
figure('Name','Maksimum IDR','Color','w');
plot(story_ids, IDR0,'k-o','LineWidth',1.4); hold on;
plot(story_ids, IDR_lin,'b-s','LineWidth',1.1);
plot(story_ids, IDR_orf,'r-d','LineWidth',1.0);
grid on; xlabel('Kat'); ylabel('Maks IDR [Delta x/h]');
legend('Dampersiz','Lineer damper','Orifisli damper','Location','best');

%% Kısa özet metni (ekrana yazdırma kaldırıldı)