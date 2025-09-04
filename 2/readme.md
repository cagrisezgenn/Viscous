# 2 Klasörü

Bu klasör, 10 katlı bir çerçevenin sismik davranışını ileri düzey damper modelleriyle inceleyen MATLAB/Octave betiklerini içerir. `1/` klasöründeki temel sürümü genişleterek, ham veya ölçeklenmiş yer hareketi kayıtlarını yükleme, kendi kendine kontrol amaçlı hesaplamalar ve ek tanılama çıktıları sağlar.

## Dosyalar

### `acc_matrix.mat`
Deprem ivmesi kayıtlarını içeren MAT dosyasıdır. `acc_matrix*` değişkenleri iki sütundan oluşur:
1. **Zaman [s]** – Deprem kaydının zaman vektörü.
2. **İvme [m/s²]** – Zamanla değişen yer ivmesi.
Bu veri seti `load_ground_motions.m` tarafından okunur ve diğer betiklere aktarılır.

### `parametreler.m`
Yapısal, damper ve akış/termal parametrelerinin tamamını tanımlar. Üç ana bölümden oluşur:
1. **Yapı parametreleri:** Kat sayısı, kütle, rijitlik, sönüm ve kat yüksekliği gibi değerlerle kütle (`M`), rijitlik (`K`) ve sönüm (`C0`) matrislerini oluşturur.
2. **Damper geometrisi ve malzeme özellikleri:** Piston çapı, orifis boyutları, viskozite gibi girdilerden seri damper sertliği (`k_sd`) ve başlangıç laminer sönüm (`c_lam0`) hesaplanır.
3. **Orifis ve termal model:** Orifis akışı için katsayılar (`orf` yapısı), termal döngü parametreleri (`thermal` yapısı), yağ/çelik özgül ısıları ve sınır değerleri tanımlanır.
Bu dosya doğrudan çalıştırılmaz; diğer betikler tarafından `parametreler;` komutu ile içeri alınır.

### `load_ground_motions.m`
`acc_matrix.mat` içindeki tüm yer hareketi kayıtlarını okur, küçük bir yüksek geçirgen filtre uygular ve kayıtları PSA(T₁, %5 sönüm) şiddetine göre ölçeklendirebilir. İşlev:
- Ham kayıtları tablo olarak döndürür.
- `T1` verildiğinde her kayıt için pseudo spektral ivme hesaplar ve hedef değere ölçekler.
- `calc_psa.m` fonksiyonunu kullanır.

### `calc_psa.m`
Tek serbestlik dereceli lineer osilatör için pseudo spektral ivmeyi hesaplar. `load_ground_motions.m` tarafından çağrılır.

### `damperlinon.m`
Ana senaryo betiğidir. Deprem ivmesi girdisini okuyarak üç senaryoyu çözer:
- **Dampersiz model** (`lin_MCK` fonksiyonu)
- **Lineer damper** (`use_orifice=false`)
- **Orifisli damper** (`use_orifice=true`), isteğe bağlı termal döngü (`use_thermal=true`)

Ek özellikler:
- `use_scaled` anahtarıyla ölçekli veya ham kayıt seçimi.
- `load_ground_motions` ile kayıtların yüklenmesi ve ölçeklenmesi.
- Birinci mod sönüm oranını hesaplayan kendi kendine kontrol bölümü.
Sonuçlar `grafik.m` betiği ile görselleştirilir.

### `grafik.m`
`damperlinon.m` tarafından çağrılan yardımcı betiktir. Üç senaryonun sonuçlarını şu şekilde sunar:
- 10. kat yer değiştirme ve mutlak ivme zaman geçmişleri
- Katlar arası maksimum göreli ötelemeler (IDR)
- Komut satırına kısa özet istatistikleri yazdırma

### `diagnostic.m`
Orifis ve termal etkiler **açık** iken damper modelini çalıştırarak her kat için tanı ölçütleri üretir. Hesaplanan büyüklükler (drift, debi, basınç kaybı, enerji vb.) `out/diagnostic.csv` dosyasına yazılır. Girdi depremi `load_ground_motions` üzerinden alınır.

## Kullanım
1. MATLAB veya Octave ortamında çalışın.
2. Ana senaryo karşılaştırmaları için `damperlinon.m` betiğini çalıştırın. `use_scaled` anahtarı ile ölçekli (`true`) veya ham (`false`) kayıt seçin.
3. Damper performansı hakkında ayrıntılı istatistikler almak için `diagnostic.m` betiğini çalıştırın. Çıktılar `out/` klasöründe `diagnostic.csv` olarak oluşturulur.

## Notlar
- Betikler, deprem verisini `acc_matrix.mat` dosyasından okur; dosya yolu betiklerle aynı klasörde olmalıdır.
- `load_ground_motions.m`, Signal Processing Toolbox yoksa yalnızca ortalama ve trend giderme işlemleri uygular.
- Termal döngü ve orifis etkilerini etkinleştirmek veya devre dışı bırakmak için `use_orifice` ve `use_thermal` anahtarları ayarlanabilir.
- Özellikler ve fonksiyonlar hakkında ayrıntılı açıklamalar için ilgili `.m` dosyalarındaki yorumlara bakın.