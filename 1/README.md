# 1 Klasörü

Bu klasör, 10 katlı bir çerçevenin sismik davranışını doğrusal ve orifisli damper modelleriyle inceleyen MATLAB betiklerini içerir. Betikler, deprem ivmesi verisini kullanarak farklı damper senaryolarını karşılaştırır ve sonuçları grafikler veya tablolara döker.

## Dosyalar

### `acc_matrix.mat`
Deprem ivmesi kayıtlarını içeren MAT dosyasıdır. `acc_matrix7` değişkeni iki sütundan oluşur:
1. **Zaman [s]**  – Deprem kaydının zaman vektörü.
2. **İvme [m/s²]** – Zamanla değişen yer ivmesi.
Bu veri seti hem `damperlinon.m` hem de `diagnostic.m` tarafından yüklenir.

### `parametreler.m`
Yapısal, damper ve akış/termal parametrelerinin tamamını tanımlar. Üç ana bölümden oluşur:
1. **Yapı parametreleri:** Kat sayısı, kütle, rijitlik, sönüm ve kat yüksekliği gibi değerlerle kütle (`M`), rijitlik (`K`) ve sönüm (`C0`) matrislerini oluşturur.
2. **Damper geometrisi ve malzeme özellikleri:** Piston çapı, orifis boyutları, viskozite gibi girdilerden seri damper sertliği (`k_sd`) ve başlangıç laminer sönüm (`c_lam0`) hesaplanır.
3. **Orifis ve termal model:** Orifis akışı için katsayılar (`orf` yapısı), termal döngü parametreleri (`thermal` yapısı), yağ/çelik özgül ısıları ve sınır değerleri tanımlanır.
Bu dosya doğrudan çalıştırılmaz; diğer betikler tarafından `parametreler;` komutu ile içeri alınır.

### `damperlinon.m`
Klasörün ana betiğidir. Deprem ivmesi girdisini okuyarak üç senaryoyu çözer:
- **Dampersiz model** (`lin_MCK` fonksiyonu)
- **Lineer damper** (`use_orifice=false`)
- **Orifisli damper** (`use_orifice=true`), isteğe bağlı termal döngü (`use_thermal=true`)

Betik önce `parametreler.m` dosyasını yükler, ardından seçilen senaryolar için `mck_with_damper` fonksiyonunu çağırır. Sonuçların özetlenmesi ve karşılaştırma grafikleri için `grafik.m` betiğini çalıştırır.

### `grafik.m`
`damperlinon.m` tarafından çağrılan yardımcı betiktir. Üç senaryonun sonuçlarını şu şekilde sunar:
- 10. kat yer değiştirme ve mutlak ivme zaman geçmişleri
- Katlar arası maksimum göreli ötelemeler (IDR)
- Komut satırına kısa özet istatistikleri yazdırma

### `diagnostic.m`
Orifis ve termal etkiler **açık** iken damper modelini çalıştırarak her kat için tanı ölçütleri üretir. Hesaplanan büyüklükler (drift, debi, basınç kaybı, enerji vb.) `out/diagnostic.csv` dosyasına yazılır. Bu betik, `damperlinon.m`'de kullanılan `mck_with_damper` fonksiyonunun tanılama amaçlı bir türetilmiş sürümünü içerir.

## Kullanım
1. MATLAB veya Octave ortamında çalışın.
2. Ana senaryo karşılaştırmaları için `damperlinon.m` betiğini çalıştırın.
3. Damper performansı hakkında ayrıntılı istatistikler almak için `diagnostic.m` betiğini çalıştırın. Çıktılar `out/` klasöründe `diagnostic.csv` olarak oluşturulur.

## Notlar
- Betikler, deprem verisini `acc_matrix.mat` dosyasından okur; dosya yolu betiklerle aynı klasörde olmalıdır.
- Termal döngü ve orifis etkilerini etkinleştirmek veya devre dışı bırakmak için `use_orifice` ve `use_thermal` anahtarları ayarlanabilir.
