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

**Sessiz Mod (quiet) ve Dışa Aktarım (do_export)**
- **quiet:** Konsol çıktıları `opts.quiet=true` iken bastırılır. `run_batch_windowed.m` ve `run_policies_windowed.m` bu bayrağı destekler.
- **do_export:** Değer `false` olduğunda dosya yazımları/`diary`/figür üretimi yapılmaz. GA değerlendirmeleri bu modu kullanır.

**Hızlı GA Sürücüsü**
- **Amaç:** Çoklu yer hareketi ve µ-ağırlıklı özet üzerinde sönümleyici tasarımını NSGA-II (`gamultiobj`) ile hızlıca optimize etmek.
- **Kaba Izgara:** Karar değişkenleri endüstriyel olarak makul adımlara “snap” edilir; `n_orf` tam sayıdır.
- **Hafif Değerlendirme:** Her değerlendirmede çizim/CSV/snapshot/diary yazımı yapılmaz; tekrar eden tasarımlar önbellekten gelir.

**Karar Vektörü ve Izgaralar**
- **x:** `[d_o_mm, n_orf, PF_tau, PF_gain, Cd0, CdInf, p_exp, Lori_mm, hA_W_perK, Dp_mm, d_w_mm, D_m_mm, n_turn, mu_ref]`
- **Adımlar:** `d_o_mm:0.1`, `n_orf:integer`, `PF_tau:0.01`, `PF_gain:0.02`, `Cd0/CdInf:0.01`, `p_exp:0.05`, `Lori_mm:1`, `hA_W_perK:25`, `Dp_mm:1`, `d_w_mm:0.5`, `D_m_mm:5`, `n_turn:integer`, `mu_ref:0.05`.

**Yeni Dosyalar**
- `quantize_step.m:1`: Bir skaler/vektörü `step*round(x/step)` ile ızgaraya oturtur.
- `decode_params_from_x.m:1`: `params` kopyalar; `x`’e göre `P.orf.d_o`, `P.n_orf`, `P.A_o`, `P.Qcap_big`, 3-bölge `P.toggle_gain` ve varsa `P.cfg.PF.(tau,gain)` günceller.
- `eval_design_fast.m:1`: GA amaç sarmalayıcı; değişkenleri kuantalize eder, `run_batch_windowed`’i `do_export=false, quiet=true` ile çağırır, µ-ağırlıklı `PFA_w` ve `IDR_w` ortalamalarını amaç yapar, QC eşiklerine yumuşak ceza ekler, sonuçları `containers.Map` ile önbellekler.
- `run_ga_driver.m:1`: NSGA-II sürücüsü; `IntCon=2`, paralel açık, çizim yok. Yalnız `out/ga_*/ga_front.(csv|mat)` yazar.

**Çalıştırma**
- **Önkoşul:** MATLAB/Global Opt. Toolbox; yol ekli: `addpath(genpath(pwd))` veya `setup`.
- **Komut:** `[X,F] = run_ga_driver(scaled, params);`
- **Çıktılar:** `out/ga_YYYYMMDD_HHMMSS/ga_front.csv` ve `ga_front.mat`.

**Kabul Kriterleri**
- **Kuantizasyon:** Belirtilen adımlara “snap”; `n_orf` tam sayı.
- **Sessizlik:** GA değerlendirmelerinde figür/CSV/diary üretilmez.
- **Önbellek:** Aynı `x` tekrar hesaplanmaz.
- **Paralel:** `gamultiobj` paralel çalışır.

## \u03bc-Robustluk (Adım 4)
Bu klasördeki `run_one_record_windowed.m` ve `run_batch_windowed.m` betikleri, s\u00f6n\u00fcmler ve yağ viskozitesindeki belirsizlikleri hesaba katan yeni bir \u03bc tarama yetene\u011fi ile g\u00fcncellendi.

- **\u03bc taramas\u0131:** `opts.mu_factors` ve `opts.mu_weights` se\u00e7enekleri ile \u03bc katsay\u0131s\u0131 \u00e7arpanlar\u0131 ve a\u011f\u0131rl\u0131klar\u0131 tan\u0131mlanabilir. Varsay\u0131lan de\u011ferler s\u0131ras\u0131yla `[0.75 1.00 1.25]` ve `[0.2 0.6 0.2]`'dir.
- **Tek pencereli \u00e7\u00f6z\u00fcm:** Arias penceresi sadece bir kez olu\u015fturulur ve t\u00fcm \u03bc senaryolar\u0131 i\u00e7in yeniden kullan\u0131l\u0131r.
- **Metrik \u00f6zetleri:** Her \u03bc senaryosu i\u00e7in pencere-i\u00e7i metrikler hesaplan\u0131r; a\u011f\u0131rl\u0131kl\u0131 ve en k\u00f6t\u00fc durum \u00f6zetleri elde edilir.
- **QC kontrolleri:** Cavitation, bas\u0131n\e7 ve s\u0131cakl\u0131k e\u015fikleri her \u03bc senaryosu i\u00e7in ayr\u0131 ayr\u0131 do\u011frulan\u0131r ve toplam sonu\u00e7 `qc_all_mu` alan\u0131nda raporlan\u0131r.
- **Toplu analiz:** `run_batch_windowed` art\u0131k nominal, a\u011f\u0131rl\u0131kl\u0131 ve en k\u00f6t\u00fc durum metriklerini i\u00e7eren geni\u015f bir \u00f6zet tablosu (`summary.table`) d\u00f6nd\u00fcr\u00fcr ve en k\u00f6t PFA/IDR de\u011ferlerinin hangi kay\u0131t ve \u03bc de\u011ferine ait oldu\u011funu g\u00fcnl\u00fc\u011fe yazar.

Bu eklemeler, damper tasar\u0131m\u0131n\u0131 viskozite sapmalar\u0131 kar\u015f\u0131s\u0131nda daha g\u00fc\u00e7l\u00fc hale getirmek i\u00e7in \u00f6ng\u00f6r\u00fclm\u00fc\u015f t\u00fcm senaryolar\u0131n birlikte de\u011ferlendirilmesini sa\u011flar.

## Is\u0131l Reset Politikalar\u0131 + QC ve K\u0131yas (Ad\u0131m 5 & 6)

Fizik/denklemler, \u00e7\u00f6z\u00fcc\u00fc ve IM/Arias mant\u0131\u011f\u0131 **de\u011fi\u015ftirilmeden**, orkestrasyon, toplama, QC ve d\u0131\u015fa aktar\u0131m katmanlar\u0131 a\u015fa\u011f\u0131daki \u015fekilde genisletildi:

- **Is\u0131l reset politikalar\u0131:** `each` | `carry` | `cooldown` (bkz. `opts.cooldown_s`, vars: 60 s).
- **Kay\u0131t s\u0131ras\u0131:** `natural` | `worst_first` | `random` (ops.). `worst_first` baz ko\u015fta belirlenen bir skora g\u00f6re s\u0131ralar.
- **\u03bc-robust a\u011f\u0131rl\u0131kl\u0131 ve en k\u00f6t\u00fc \u00f6zetler:** \u03bc senaryolar\u0131 \u00fczerinden a\u011f\u0131rl\u0131kl\u0131 ortalama ve metrik-uygun en k\u00f6t durum (\u00f6rn. `mu_end` i\u00e7in min) hesaplan\u0131r.
- **QC e\u015fikleri (konfig\u00fcre edilebilir):** `opts.thr = struct('dP95_max',50e6,'Qcap95_max',0.5,'cav_pct_max',0,'T_end_max',75,'mu_end_min',0.5)`.
- **Politika k\u0131yas kural\u0131:** Baz `each/natural` ko\u015funa g\u00f6re `|\u0394PFA_w|, |\u0394IDR_w| \le 0.15 \times` baz ortalama.
- **Zengin log + eksport:** `out/<timestamp>/` alt\u0131na `console.log`, endeks/\u00f6zet/policy CSV'leri ve geni\u015f `snapshot.mat` yaz\u0131l\u0131r.

### Yeni/Genel Se\u00e7enekler (opts)

- `opts.policies`: `{'each','carry','cooldown'}` alt k\u00fcmesi.
- `opts.orders`: `{'natural','worst_first','random'}` alt k\u00fcmesi.
- `opts.cooldown_s_list`: `cooldown` test s\u00fcreleri (\u00f6rn. `[60 180 300]`).
- `opts.rng_seed`: `random` s\u0131ralama i\u00e7in.
- `opts.rank_metric`: `E_orifice_win` (vars) | `pfa_top` | `idr_max` | `pfa_w`.
- `opts.mu_factors`, `opts.mu_weights`: \u03bc taramas\u0131.
- `opts.thr`: QC e\u015fikleri.
- `opts.do_export`: `run_policies_windowed` i\u00e7in varsay\u0131lan `true`.

### Dosyalar ve Roller

- `run_one_record_windowed.m`: Tek kay\u0131t ko\u015fumlar\u0131. Yeni alanlar: `T_start`, `T_end`, `mu_end`, `clamp_hits`, `qc_all_mu`, `PFA_top`, `IDR_max`, `dP_orf_q95`, `Qcap_ratio_q95`, `cav_pct`, `t5`, `t95`, `coverage`. A\u011f\u0131rl\u0131kl\u0131/en k\u00f6t\u00fc \u00f6zetlerde `cav_pct` ve `mu_end` de kapsan\u0131r.
- `run_batch_windowed.m`: Toplu ko\u015fum ve \u00f6zet tablo. Eklenen s\u00fctunlar: `policy`, `order`, `cooldown_s`, `PFA_w`, `IDR_w`, `PFA_worst`, `IDR_worst`, `dP95_worst`, `Qcap95_worst`, `cav_pct_worst`, `T_end_worst`, `mu_end_worst`, `qc_all_mu`, `T_start`, `T_end`, `mu_end`, `clamp_hits`, `t5`, `t95`, `coverage`, `rank_score` (`worst_first` d\u0131\u015f\u0131nda `NaN`).
- `run_policies_windowed.m`: Policy/order/(cooldown) kombinasyonlar\u0131n\u0131 \u00e7al\u0131\u015ft\u0131r\u0131r, baz (each/natural) ile kar\u015f\u0131la\u015ft\u0131rma yapar, tek sat\u0131rl\u0131 \u00f6zetler basar ve \u00e7\u0131kt\u0131lar\u0131 d\u0131\u015fa aktar\u0131r.
- `export_results.m`: D\u0131\u015fa aktar\u0131mlar. `scaled_index.csv` (ek kolonlar: `s_clipped`, `trimmed`), `summary.csv`, `policy_index.csv`, her kombinasyona `policy_*.csv`, ve geni\u015f `snapshot.mat` (`IM_mode`, `band_fac`, `s_bounds`, `TRIM_names`, ve m\u00fcsaitse `params_derived.C_th`).

### Konsol/Log (console.log)

- Ko\u015f ba\u015fl\u0131\u011f\u0131: zaman damgas\u0131, `outdir`, IM ayarlar\u0131, politikalar/s\u0131p, `cooldown_s_list`, `rng_seed`, varsa `TRIM` listesi, QC e\u015fikleri.
- Baz sat\u0131r\u0131: `PFA_w`, `IDR_w` ort., `dP95_worst` (MPa), `Qcap95_worst`, `cav%_worst`, `T_end_worst` (C), `mu_end_worst` (min), `qc_rate=pass/total`.
- `worst_first` s\u0131ralamas\u0131: `opts.rank_metric` ve s\u0131ral\u0131 kay\u0131t adlar\u0131.
- Her kombinasyon i\u00e7in: `policy`, `order`, `cd`, `PFA_w` ve `IDR_w` (baz ortalamaya g\u00f6re %\u0394), `T_end_worst` (maks), `qc_rate=pass/total`, `%15` kural\u0131 i\u00e7in `OK/FAIL`. Opsiyonel: `clamp_hits` toplam\u0131 ve isim listesi.

### D\u0131\u015fa Aktar\u0131mlar (out/<timestamp>/)

- `console.log`: T\u00fcm \u00e7al\u0131\u015fma g\u00fcnl\u00fc\u011f\u00fc.
- `scaled_index.csv`: `name, dt, dur, PGA, PGV, IM, scale, s_clipped, trimmed`.
- `summary.csv`: Baseline `summary.table` (yukar\u0131daki s\u00fctunlar).
- `policy_index.csv`: `policy, order, cooldown_s, qc_rate, PFA_w_mean, IDR_w_mean, dP95_worst_max, T_end_worst_max`.
- `policy_<policy>_<order>_cd<sec>.csv`: Her kombinasyon i\u00e7in tam tablo.
- `snapshot.mat`: `params, opts, scaled, P, IM_mode, band_fac, s_bounds, TRIM_names, params_derived`.

### H\u0131zl\u0131 Ba\u015flang\u0131\u00e7

- \u00d6nerilen ayarlar:
  - `opts.policies = {'each','carry','cooldown'}`
  - `opts.orders = {'natural','worst_first'}`
  - `opts.cooldown_s_list = [60 180 300]`
  - `opts.mu_factors = [0.75 1.00 1.25]`, `opts.mu_weights = [0.2 0.6 0.2]`
  - `opts.rank_metric = 'E_orifice_win'`
  - `opts.thr = struct('dP95_max',50e6,'Qcap95_max',0.5,'cav_pct_max',0,'T_end_max',75,'mu_end_min',0.5)`
  - \u00c7al\u0131\u015ft\u0131r: `P = run_policies_windowed(scaled, params, opts);`

### Ek Notlar

- Fizik ve IM/Arias mant\u0131\u011f\u0131 de\u011fi\u015fmedi; t\u00fcm yenilikler orkestrasyon ve IO katmanlar\u0131ndad\u0131r.
- `rank_score` yaln\u0131zca `order='worst_first'` i\u00e7in anlaml\u0131 olup di\u011ferlerinde `NaN` yaz\u0131l\u0131r.
- `policy_index.csv`, kombinasyon \u00e7\u0131kt\u0131lar\u0131ndan (mean/maks) t\u00fcretilmi\u015f h\u0131zl\u0131 k\u00fcme \u00f6zetleri sa\u011flar.
