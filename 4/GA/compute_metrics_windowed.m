function metr = compute_metrics_windowed(t, x, a_rel, ag, ts, story_height, win, params)
%COMPUTE_METRICS_WINDOWED Zaman penceresi içindeki tepkileri hesaplar.
%   METR = COMPUTE_METRICS_WINDOWED(T,X,A_REL,AG,TS,STORY_HEIGHT,WIN,PARAMS)
%   fonksiyonu, WIN.IDX tarafından tanımlanan zaman aralığında yapının
%   performansına ilişkin metrikleri üretir. METR değişkeni tepe kat mutlak
%   ivmesi, katlar arası ötelenme oranları, orifis basınç istatistikleri,
%   enerji ölçümleri ve nihai ("sıcak") damper katsayısına bağlı modal sönüm
%   oranlarını içerir.
%
%   Girdi değişkenleri T, X, A_REL ve AG sırasıyla zaman vektörü, kat yer
%   değiştirmeleri, göreli kat ivmeleri ve yer ivmesini temsil eder. TS
%   yapısal analizin ürettiği ek zaman serilerini içerir; bu yapı sadece
%   metrik hesapları için gerekli dizi alanlarını barındırır. STORY_HEIGHT her
%   katın yüksekliğidir. WIN.IDX ilgilenilen pencereyi seçen mantıksal
%   vektördür. PARAMS yapısal ve damper parametrelerini içerir.

idx = win.idx;

%% Temel Tepki
% Bu bölüm, tepe kat ivmesi ve katlar arası göreli hareketler gibi yapısal
% yanıtın temel göstergelerini hesaplar.
% Tepe katın mutlak ivmesi
a_top_abs = a_rel(:,end) + ag(:);
metr.PFA = max(abs(a_top_abs(idx)));

% Pencere içindeki tepe kat mutlak yerdeğiştirmesi ve ivmesi
metr.x10_max_D    = max(abs(x(idx,end)));
metr.a10abs_max_D = max(abs(a_top_abs(idx)));

% Katlar arası göreli ötelenme oranının maksimumu
drift = (x(:,2:end) - x(:,1:end-1)) / story_height;
metr.IDR = max(max(abs(drift(idx,:))));

%% Kat Bazlı İstatistikler
% Her kat için yüzde 50 ve yüzde 95'lik değerler ile kavitasyon oranı
% gibi istatistikler hesaplanır.

abs_dP = abs(ts.dP_orf(idx,:));
abs_Q  = abs(ts.Q(idx,:));
Qcap_ratio = abs_Q ./ max(params.Qcap_big, eps);
abs_story_force = abs(ts.story_force(idx,:));

% Reynolds sayısının pencere içindeki maksimumu
try
    if all(isfield(params,{'Qcap_big','Ap','Ao','rho','orf'})) && ...
            isfield(ts,'mu') && isfield(ts,'dvel') && ...
            all(isfield(params.orf,{'veps','d_o'}))
        dvel_win = ts.dvel(idx,:);
        qmag = params.Qcap_big * tanh((params.Ap/params.Qcap_big) * ...
            sqrt(dvel_win.^2 + params.orf.veps^2));
        mu_win = ts.mu(idx);
        mu_mat = repmat(mu_win,1,size(qmag,2));
        Ao = params.Ao; if numel(Ao)==1, Ao = Ao*ones(1,size(qmag,2)); end
        Ao_mat = repmat(Ao,size(qmag,1),1);
        Re = (params.rho .* qmag ./ max(Ao_mat .* mu_mat,1e-9)) .* ...
             max(params.orf.d_o,1e-9);
        metr.Re_max = max(Re,[],'all');
    else
        metr.Re_max = NaN;
    end
catch
    metr.Re_max = NaN;
end

% 50. ve 95. yüzdelik değerler
dP_q50         = local_quantile(abs_dP, 0.50);
dP_q95         = local_quantile(abs_dP, 0.95);
Q_q50          = local_quantile(abs_Q, 0.50);
Q_q95          = local_quantile(abs_Q, 0.95);
Qcap95 = local_quantile(Qcap_ratio, 0.95);
story_force_q95= local_quantile(abs_story_force, 0.95);

% Her kat için ortalama kavitasyon yüzdesi
cav_mean = mean(ts.cav_mask(idx,:),1);

% Hikaye kesme kuvvetine göre kritik katın belirlenmesi
[metr.story_force_q95, metr.which_story] = max(story_force_q95);
ws = metr.which_story;

% Kritik katın istatistiklerinin saklanması
metr.dP95   = dP_q95(ws);
metr.dP50   = dP_q50(ws);
metr.Q_q95  = Q_q95(ws);
metr.Q_q50  = Q_q50(ws);
metr.Qcap95 = Qcap95(ws);
metr.cav_pct         = cav_mean(ws);

%% Enerji Hesapları
% Pencerede ve tüm süreçte biriken enerji bileşenleri değerlendirilir.
w_first = find(idx,1,'first');
w_last  = find(idx,1,'last');
i0 = max(w_first-1,1);

% Toplam süreç sonundaki orifis ve yapı enerjileri
metr.E_orifice = ts.E_orf(end);
metr.E_struct  = ts.E_struct(end);
metr.E_ratio   = metr.E_orifice / max(metr.E_struct, eps);
% Toplam enerji ve ortalama mekanik güç
metr.energy_tot = metr.E_orifice + metr.E_struct;
try
    if isfield(ts,'P_sum') && ~isempty(ts.P_sum)
        metr.P_mech = mean(ts.P_sum(idx));
    else
        metr.P_mech = NaN;
    end
catch
    metr.P_mech = NaN;
end

% Seçilen pencere içindeki enerji birikimleri
metr.E_orifice_win = ts.E_orf(w_last) - ts.E_orf(i0);
metr.E_struct_win  = ts.E_struct(w_last) - ts.E_struct(i0);
metr.E_ratio_win   = metr.E_orifice_win / max(metr.E_struct_win, eps);

%% Termal Metrikler
% Yağ viskozitesi gibi termal büyüklükler değerlendirilir.
if isfield(ts,'mu')
    metr.mu_end = ts.mu(w_last);
else
    metr.mu_end = NaN;
end

% ---------------- Sıcak viskozite ile modal sönüm -------------------
req_fields = {'M','K','C0','k_sd'};
if all(isfield(params,req_fields)) && isfield(ts,'c_lam')
    M  = params.M;  K = params.K;  C0 = params.C0;
    k_sd = params.k_sd;
    nStories = size(M,1) - 1;
    c_lam = ts.c_lam;

    Kadd = zeros(size(M));
    Cadd = zeros(size(M));
    for i=1:nStories
        idx2 = [i i+1];
        % R katsayısı lineer olarak kullanılır (R), karesi (R^2) alınmaz
        k_eq = k_sd;
        c_eq = c_lam;
        kM = k_eq * [1 -1; -1 1];
        cM = c_eq * [1 -1; -1 1];
        Kadd(idx2,idx2) = Kadd(idx2,idx2) + kM;
        Cadd(idx2,idx2) = Cadd(idx2,idx2) + cM;
    end
    Ktot = K + Kadd;
    Ctot = C0 + Cadd;
    [V,D] = eig(Ktot,M);
    [w2,ord] = sort(diag(D),'ascend');
    V = V(:,ord);
    phi1 = V(:,1); phi2 = V(:,2);
    w1 = sqrt(w2(1)); w2s = sqrt(w2(2));
    n1 = phi1.'*M*phi1; n2 = phi2.'*M*phi2;
    z1 = (phi1.'*Ctot*phi1)/(2*w1*n1);
    z2 = (phi2.'*Ctot*phi2)/(2*w2s*n2);
    metr.zeta1_hot = z1;
    metr.z2_over_z1_hot = z2 / max(z1, eps);
else
    metr.zeta1_hot = NaN;
    metr.z2_over_z1_hot = NaN;
end

% ----------------------- PF metrikleri (isteğe bağlı) -----------------
try
    if isfield(ts,'PF')
        PF_abs = abs(ts.PF(idx,:));
    metr.PF_p95 = local_quantile(PF_abs(:,ws), 0.95);
    end
catch
end

% Parametrelerden gerekli alanlar kayıt altına alınır
param_fields = {'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref'};
for ii = 1:numel(param_fields)
    fn = param_fields{ii};
    if isfield(params, fn)
        metr.(fn) = params.(fn);
    end
end

end

function q = local_quantile(A, p)
%LOCAL_QUANTILE Verilen verinin p yüzdelik değerini hesaplar.
%   A matrisi ve 0-1 aralığındaki p değeri ile bu yardımcı fonksiyon
%   MATLAB'in QUANTILE fonksiyonunu çağırarak yüzdelik değeri döndürür.
q = quantile(A, p);
end

