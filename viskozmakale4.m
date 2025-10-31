%% GA Sürücüsü ve Kurulum
% RUN_GA_DRIVER hibrit genetik algoritma (GA) iş akışını tek dosyada yönetir.
% Amaç ve Kapsam: Bu ana sürücü, tasarım değişkenlerinin sınırları içindeki
% kantitatif arama problemini çözer, hızlı ceza hesaplarını tetikler ve her
% Pareto çözümü için pencere tabanlı metrikleri toplar.
% Girdiler:
%   scaled [-]: Önceden ölçeklendirilmiş yer hareketi kayıtları yapısı.
%   params [-]: Başlangıç damper ve yapı parametreleri (SI birimleri).
%   optsEval [-]: Ceza çekirdeği, QC eşiği ve pencereleme ayarları.
%   optsGA [-]: Popülasyon, mutasyon ve paralel ayarlarını içeren GA opsiyonları.
% Çıktılar:
%   X [-]: Pareto ön cephesindeki karar vektörü satırları (mm değerleri henüz
%        kuantize edilmemiş ham formdadır).
%   F [-]: [f_pen, f1, f2] amaç değerlerini içeren tablo; f_pen yumuşak ceza
%        terimidir, f1 ve f2 ise performans göstergeleridir.
%   gaout [-]: GA yürütme meta verisi (exitflag, iterasyon istatistikleri).
% Varsayımlar ve Sınırlar:
%   - Parametreler `prepare_inputs` ile üretildiğinde tutarlı birimler sağlanır.
%   - Parpool açılımı başarısız olursa seri hesap kabul edilir.
%   - Karar vektörü sınırları mm bazındadır; değerlendirme sırasında mm→m
%     dönüşümü alt fonksiyonlarda yapılır.
% Ölçüler/Birimler: Q [m^3/s], Δp [Pa], P_mech_sum [J], energy_tot_sum [J],
%   T_end [°C], mu_end [Pa·s], cav_pct [%], IDR [-]:, PFA [-]:.
% Yöntem Özeti: GA sürücüsü Arias yoğunluğu penceresi (t5–t95) ile
% pencereleme yapan `run_batch_windowed` fonksiyonunu çağırır. Termal durum
% modları (her kayıt sonrası reset, taşıma veya soğuma) `run_one_record_windowed`
% içindeki politika yorumlarında açıklanır. Güç-enerji ilişkisi P_mech ≈
% Σ(Δp·q·Δt) olup [W·s] ≙ [J]; bu integral `summarize_metrics_table`
% tarafından tutulur. Yumuşak ceza `pen_raw = Σ W_i·δ_i` ve `pen = λ·pen_raw`
% şeklinde uygulanır; δ_i göreli sapma, λ ölçek katsayısıdır ve hepsi boyutsuzdur.
function [X,F,gaout] = run_ga_driver(scaled, params, optsEval, optsGA)

narginchk(0,4);

if nargin < 1 || isempty(scaled),  scaled  = []; end
if nargin < 2 || isempty(params),  params  = []; end
if nargin < 3 || isempty(optsEval), optsEval = struct; end
if nargin < 4 || isempty(optsGA),   optsGA   = struct; end

% Gerekli girdiler sağlanmadıysa otomatik hazırla.
if isempty(scaled) || isempty(params)
    [scaled, params] = prepare_inputs(optsGA);
end
assignin('base','scaled',scaled);
assignin('base','params',params);
% === Parpool açılışı (temizlik + iş parçacığı sınırı) ===
usePool = true;
try
    usePool = parpool_hard_reset(16);
catch ME
    warning('run_ga_driver:parpool', 'Parallel pool unavailable: %s', ME.message);
    usePool = false;
end

assert(~isempty(scaled), 'run_ga_driver: scaled dataset is empty.');
assert(~isempty(params), 'run_ga_driver: params is empty.');

% QC eşikleri thr [Pa, -, °C, Pa·s] tek kaynaktan toplanır ve parametre
% bağımlı eksikler compute_thr ile tamamlanır; kullanıcı thr alanı baskın kalır.
thr_src = resolve_thr_base(optsEval);
thr = compute_thr(params, thr_src);
optsEval.thr = thr;

%% Amaç Fonksiyonu ve Ceza Çekirdeği
% Penalty çekirdeği [boyutsuz] parametreleri λ [-]:, kuvvet pwr [-]: ve ağırlık
% W_i [-]: hidrolik kritikler (Δp, Qcap, cav, T, μ) için burada set edilir. pen_raw =
% Σ W_i·δ_i ifadesi ile toplanır ve pen = λ·pen_raw olarak ölçeklenir; δ_i
% göreli sapmadır.
%% GA Sürücüsü ve Kurulum — Optimizasyon Ayarları
% GA amaç fonksiyonu ve optimizasyon seçeneklerini hazırla.
    rng(42);

    % Karar vektoru: [d_o_mm, n_orf, Cd0, CdInf, p_exp, Lori_mm, hA_W_perK, Dp_mm, d_w_mm, D_m_mm, n_turn, mu_ref]
    lb = [1.0, 4, 0.60, 0.75, 0.90, 120, 200, 120, 10,  90,  6, 0.80];
    ub = [3.0, 8, 0.90, 1.00, 1.50, 200, 600, 240, 16, 160, 18, 1.5];
    IntCon = [2 11];  % n_orf ve n_turn tam sayı

    obj = @(x) eval_design_fast(x, scaled, params, optsEval); % içerde kuantize/clamplar

  % --- HIZLI ÖN TARAMA AYARLARI ---
POP_SIZE_QUICK = 40;   % Düşük popülasyon yeterli
MAX_GEN_QUICK = 1;    % Çok az jenerasyon (Belki 15-20 de olabilir)

options = optimoptions('gamultiobj', ...
       'PopulationSize',    POP_SIZE_QUICK, ... % Düşük
       'MaxGenerations',    MAX_GEN_QUICK, ...  % Düşük
       'CrossoverFraction', 0.8, ...
       'MutationFcn',       {@mutationadaptfeasible}, ...
       'ParetoFraction',    0.4, ... % Çeşitliliği korumak önemli
       'StallGenLimit',     10, ... % Düşük
       'DistanceMeasureFcn','distancecrowding', ...
       'OutputFcn',         @(options,state,flag) ga_out_best_pen(options,state,flag, scaled, params, optsEval), ...
       'UseParallel',       usePool, ...
       'Display','iter', ...
       'FunctionTolerance', 1e-4); % Toleransı biraz gevşetebiliriz

        [X,F,exitflag,output] = gamultiobj(obj, numel(lb), [],[],[],[], lb, ub, [], IntCon, options);
    gaout = struct('exitflag',exitflag,'output',output);
    assignin('base','X',X);
assignin('base','F',F);
assignin('base','gaout',gaout);

    %% Sonuçların Paketlenmesi
    % GA tamamlandıktan sonra sonuçları dosyalara kaydet.
    Opost = struct('thr',thr,'penalty',penopts);
    tstamp = datestr(now,'yyyymmdd_HHMMSS_FFF');
    outdir = fullfile('out', ['ga_' tstamp]);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    date_str = tstamp;
    front = struct('X',X,'F',F,'options',options,'date_str',date_str, ga_out_best_pen);
    save(fullfile(outdir,'ga_front.mat'), '-struct', 'front', '-v7.3');

    % === Re-evaluate Pareto designs to collect metrics per-row ===
    nF = size(X,1);
    x10_max_damperli  = zeros(nF,1);  a10abs_max_damperli = zeros(nF,1);
    dP95   = zeros(nF,1);  Qcap95 = zeros(nF,1); cav_pct = zeros(nF,1);
    T_end  = zeros(nF,1);  mu_end  = zeros(nF,1);
    Q_q50 = zeros(nF,1);  Q_q95 = zeros(nF,1);
    dP50   = zeros(nF,1);
    energy_tot_sum   = zeros(nF,1);  E_orifice_sum = zeros(nF,1);   E_struct_sum  = zeros(nF,1); E_ratio = zeros(nF,1); P_mech_sum = zeros(nF,1);
    PFA_mean   = zeros(nF,1);  IDR_mean  = zeros(nF,1);

    % ceza bileşenleri (eval ile aynı)
    pen     = zeros(nF,1);
    pen_dP  = zeros(nF,1); pen_Qcap = zeros(nF,1); pen_cav = zeros(nF,1); pen_T = zeros(nF,1); pen_mu = zeros(nF,1);
    % parfor güvenliği için tüm yardımcılar subfunction'dır; nested kullanılmaz.
    parfor i = 1:nF
        Xi = quant_clamp_x(X(i,:));
        Pi = decode_params_from_x(params, Xi);
        Si = run_batch_windowed(scaled, Pi, Opost);
        metrics = summarize_metrics_table(Si.table, Opost, lambda, pwr, W);

        PFA_mean(i) = metrics.PFA_mean;
        IDR_mean(i) = metrics.IDR_mean;
        x10_max_damperli(i) = metrics.x10_max_damperli;
        a10abs_max_damperli(i) = metrics.a10abs_max_damperli;
        dP95(i)   = metrics.dP95;
        Qcap95(i) = metrics.Qcap95;
        cav_pct(i) = metrics.cav_pct;
        T_end(i)  = metrics.T_end;
        mu_end(i) = metrics.mu_end;
        Q_q50(i)  = metrics.Q_q50;
        Q_q95(i)  = metrics.Q_q95;
        dP50(i)   = metrics.dP50;
        energy_tot_sum(i) = metrics.energy_tot_sum;
        E_orifice_sum(i)  = metrics.E_orifice_sum;
        E_struct_sum(i)   = metrics.E_struct_sum;
        E_ratio(i)        = metrics.E_ratio;
        P_mech_sum(i)     = metrics.P_mech_sum;
        pen_dP(i) = metrics.pen_dP;
        pen_Qcap(i) = metrics.pen_Qcap;
        pen_cav(i) = metrics.pen_cav;
        pen_T(i)   = metrics.pen_T;
        pen_mu(i)  = metrics.pen_mu;
        pen(i)     = metrics.pen;
    end

    % Satır başına dizilerden T tablosunu oluştur
    data = [X F PFA_mean IDR_mean pen pen_dP pen_Qcap pen_cav pen_T pen_mu ...
            x10_max_damperli a10abs_max_damperli dP95 Qcap95 cav_pct T_end mu_end ...
            Q_q50 Q_q95 dP50 energy_tot_sum E_orifice_sum E_struct_sum E_ratio P_mech_sum];
    T = array2table(data, 'VariableNames', ...
       {'d_o_mm','n_orf','Cd0','CdInf','p_exp','Lori_mm','hA_W_perK','Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref', ...
        'f_pen','f1','f2','PFA_mean','IDR_mean','pen','pen_dP','pen_Qcap','pen_cav','pen_T','pen_mu', ...
        'x10_max_damperli','a10abs_max_damperli','dP95','Qcap95','cav_pct','T_end','mu_end', ...
        'Q_q50','Q_q95','dP50','energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio','P_mech_sum'});

    % === BASELINE (pre-GA) ROW: params başlangıcıyla tek koşu, ilk satır ===
    % Karşılaştırılabilirlik için baseline tasarım GA ile aynı QC [Pa, -, °C, Pa·s]
    % ve ceza çekirdeği (λ, pwr, W) altında yeniden değerlendirilir; böylece
    % CSV/MAT dosyalarında yer alan ilk satır bilimsel açıdan tutarlı kalır.
    T = prepend_baseline_row(T, params, scaled, Opost, lambda, pwr, W);

    write_pareto_results(T, outdir);
    raporla_ga_sonuclari(X, F, scaled, params, optsEval, gaout);
    ciz_ga_grafikleri(X, scaled, params, optsEval);
    % En iyi K tasarımın parametre listesi
    K = min(10, size(X,1));
    if K > 0
        idx = unique(round(linspace(1,size(X,1),K)));
        params_list = cell(numel(idx),1);
        for i = 1:numel(idx)
            params_list{i} = decode_params_from_x(params, X(idx(i),:)); %#ok<AGROW>
        end
        tmp_struct = struct('params_list',{params_list});
        save(fullfile(outdir,'ga_front.mat'), '-struct', 'tmp_struct', '-append');
    end

end

function [f, meta, details] = eval_design_fast(x, scaled, params_base, optsEval)
% EVAL_DESIGN_FAST hızlı amaç fonksiyonu değerlendirmesini yürütür.
% Amaç ve Kapsam: Kuantize edilmiş karar vektörünü (mm seviyesinde) dinamik
%   simülasyona aktarır, pencere tabanlı metrikleri toplar ve [f_pen, f1, f2]
%   vektörünü döndürür.
% Girdiler:
%   x [-]: Tasarım karar vektörü; mm temelli bileşenler quant_clamp_x ile gridlenir.
%   scaled [-]: Ölçekli yer hareketi kayıtları seti.
%   params_base [-]: Referans parametre yapısı (SI birimleri).
%   optsEval [-]: Ceza, QC, pencere ve termal politikaları içeren yapı.
% Çıktılar:
%   f [-]: [f_pen, f1, f2] amaç değerleri; f1 = mean(PFA) [m/s²], f2 = mean(IDR) [-]:.
%   meta [-]: Ceza parçaları, enerji toplamları ve QC sonuçlarını içeren özet.
%   details [-]: run_batch_windowed çıktısının isteğe bağlı tam kopyası.
% Varsayımlar ve Sınırlar:
%   - x vektörü 12 bileşenli olup orifis, termal ve malzeme parametrelerini içerir.
%   - optsEval.penalty alanı λ [-]:, pwr [-]:, W_i [-]: ve cav_free [-] değerlerini sağlar.
%   - scaled dataset_signature ile cache anahtarı oluşturulur; veri değişirse cache temizlenmelidir.
% Ölçüler/Birimler: dP95 [Pa], Qcap95 [-]:, cav_pct [%], T_end [°C], μ_end [Pa·s],
%   P_mech_sum [J], energy_tot_sum [J], Re [-]:, PFA [m/s²], IDR [-]:.
% Yöntem Özeti: Karar vektörü önce quant_clamp_x ile belirtilen ızgara
%   adımlarına (bkz. Tasarım Değişkenlerinin Izgaraya Kilitlenmesi) zorlanır.
%   decode_params_from_x mm→m dönüşümlerini uygulayarak fiziksel modele aktarır.
%   run_batch_windowed Arias yoğunluğu penceresi (t5–t95) üzerinden her kaydı
%   çözer ve pen_raw = Σ W_i·δ_i hesaplanır; pen = λ·pen_raw^(pwr) şeklinde
%   yumuşak ceza oluşturulur. cav_free [-] ölü-bölgesi cav_pct < cav_free
%   olduğunda δ_cav = 0 kabul eder; bu nedenle hafif kavitasyon ceza üretmez.
    details = struct();
    % Kuantizasyon: eval ve GA sonrası aynı ızgara [mm, -, -] (bkz. quant_clamp_x)
    x = x(:)';
    x = quant_clamp_x(x);

    

    persistent memo;
    if isempty(memo), memo = containers.Map(); end
    key = jsonencode([x, dataset_signature(scaled)]);

    O = struct();
    if nargin >= 4 && ~isempty(optsEval), O = optsEval; end

    % QC thr [Pa, -, °C, Pa·s] ve ceza çekirdeği parametreleri (λ [-]:, güç [-]:, W_i [-])
    % tek kaynaktan okunur; compute_thr ile eksikler tamamlanıp alt yordamlarla paylaşılır.
    thr_src = resolve_thr_base(O);
    thr = compute_thr(params_base, thr_src);
    O.thr = thr;

    penopt = util_getfield_default(O,'penalty', struct());
    lambda = util_getfield_default(penopt,'lambda',5);
    pwr    = util_getfield_default(penopt,'power',1.2);
    W      = util_getfield_default(penopt,'W', struct());
    if ~isstruct(W), W = struct(); end
    W.dP   = util_getfield_default(W,'dP',2);
    W.Qcap = util_getfield_default(W,'Qcap',3);
    W.cav  = util_getfield_default(W,'cav',0.25);
    W.T    = util_getfield_default(W,'T',0.5);
    W.mu   = util_getfield_default(W,'mu',0.5);
    cav_free = util_getfield_default(penopt,'cav_free',0.002);  % [-] ≈ %0.2 serbest bölge
    penopt.lambda = lambda;
    penopt.power  = pwr;
    penopt.W      = W;
    penopt.cav_free = cav_free;
    O.penalty     = penopt;

    if isKey(memo, key)
        meta = memo(key);
        f = meta.f;
        if nargout > 2
            P_cached = decode_params_from_x(params_base, meta.x);
            details = run_batch_windowed(scaled, P_cached, O);
        end
        return;
    end

    P = decode_params_from_x(params_base, x);

    % Güvenli değerlendirme (GA sırasında IO yok)
    S = run_batch_windowed(scaled, P, O);

    % --- HARD-KILL (ham Qcap95 [-]:, cav_pct [-]: ve süre eşiği [-]) ---
    qcapv_raw = enforce_finite(S.table.Qcap95);
    cavv_raw  = enforce_finite(S.table.cav_pct);

    kill_reason = '';
    if any(qcapv_raw > 0.90)
        kill_reason = 'Qcap95>0.90';
    end
    if any(cavv_raw > 0.01)
        if ~isempty(kill_reason), kill_reason = [kill_reason '+']; end
        kill_reason = [kill_reason 'cav>0.01'];
    end

    if ~isempty(kill_reason)
        f = [1e6, 1e6, 1e6];
        meta = struct('x',x,'f',f,'hard_kill',true, ...
                      'hard_kill_reason', kill_reason, ...
                      'pen',0,'pen_raw',0,'pen_dP',0,'pen_Qcap',0,'pen_cav',0,'pen_T',0,'pen_mu',0, ...
                      'qcap95_raw_max',max(qcapv_raw), 'cav_pct_max',max(cavv_raw), ...
                      'f_pen',f(1),'f1',f(2),'f2',f(3));
        fprintf('[HARD-KILL] %s\n', kill_reason);
        memo(key) = meta;
        if nargout > 2, details = S; end
        return;
    end

    f1 = mean(S.table.PFA);
    f2 = mean(S.table.IDR);
    dP95v   = enforce_finite(S.table.dP95);
    T_endv   = enforce_finite(S.table.T_end);
    mu_endv  = enforce_finite(S.table.mu_end);
    % --- PENALTY (λ [-]:, güç [-]:, W_i [-] ile boyutsuz çarpan) ---
    rel = @(v,lim) max(0,(v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0,(lim - v)./max(lim,eps)).^pwr;
    pen_dP   = isfield(thr,'dP95_max')   * mean(rel(dP95v,  thr.dP95_max));
    pen_Qcap = isfield(thr,'Qcap95_max') * mean(rel(qcapv_raw, thr.Qcap95_max));
    cav_lim = util_getfield_default(thr,'cav_pct_max',0);
    cav_eff = max(0, cavv_raw - cav_free);
    if cav_lim <= 0
        pen_cav = mean(cav_eff.^pwr);
    else
        cav_lim_eff = max(cav_lim - cav_free, eps);
        pen_cav = mean((cav_eff ./ max(cav_lim_eff, eps)).^pwr);
    end
    pen_T    = isfield(thr,'T_end_max')  * mean(rel(T_endv,  thr.T_end_max));
    pen_mu   = isfield(thr,'mu_end_min') * mean(rev(mu_endv, thr.mu_end_min));
    pen_core = W.dP*pen_dP + W.Qcap*pen_Qcap + W.cav*pen_cav + W.T*pen_T + W.mu*pen_mu;
    if ~isfinite(pen_core)
        pen_core = 0;
    end
    pen_total = lambda * pen_core;
    if ~isfinite(pen_total)
        pen_total = 1e6;
    end
    pen_total = max(pen_total, 0);
    % Ceza boyutsuzdur; lexicographic GA için f = [f_pen, f1, f2].
    f = [pen_total, f1, f2];
    pen_parts = struct('dP',pen_dP,'Qcap',pen_Qcap,'cav',pen_cav,'T',pen_T,'mu',pen_mu);
    pf = {'dP','Qcap','cav','T','mu'};
    for ii=1:numel(pf)
        fn = pf{ii};
        if ~isfield(pen_parts,fn) || ~isfinite(pen_parts.(fn)), pen_parts.(fn) = 0; end
    end
    x10_max_damperli_local = 0; a10abs_max_damperli_local = 0;
    if isfield(S,'table') && istable(S.table)
        if ismember('x10_max_damperli', S.table.Properties.VariableNames)
            x10_max_damperli_local = max(S.table.x10_max_damperli);
        end
        if ismember('a10abs_max_damperli', S.table.Properties.VariableNames)
            a10abs_max_damperli_local = max(S.table.a10abs_max_damperli);
        end
    end
    meta = struct('x',x,'f',f,'PFA_mean',f1,'IDR_mean',f2, ...
                    'pen',pen_total,'pen_raw',pen_core, ...
                    'pen_dP',pen_parts.dP,'pen_Qcap',pen_parts.Qcap, ...
                    'pen_cav',pen_parts.cav,'pen_T',pen_parts.T,'pen_mu',pen_parts.mu, ...
                   'pen_parts',pen_parts,'x10_max_damperli',x10_max_damperli_local, ...
                   'a10abs_max_damperli',a10abs_max_damperli_local, ...
                   'hard_kill',false,'hard_kill_reason','', ...
                   'qcap95_raw_max',max(qcapv_raw), 'cav_pct_max',max(cavv_raw), ...
                   'f_pen',f(1),'f1',f(2),'f2',f(3));
    % === Penaltı sürücüleri ve diagnostikleri ekle (varsa) ===
        if isfield(S,'table') && istable(S.table)
            candCols = {'dP95','Qcap95','cav_pct','T_end','mu_end', ...
                         'Q_q50','Q_q95','dP50','energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio','P_mech_sum', ...
                         'x10_max_damperli','a10abs_max_damperli','Re_max'};
            % Gerekirse alan adları için takma adlar sağla
                if ismember('E_orifice_sum', S.table.Properties.VariableNames) && ...
                   ismember('E_struct_sum', S.table.Properties.VariableNames) && ...
                   ~ismember('energy_tot_sum', S.table.Properties.VariableNames)
                    S.table.energy_tot_sum = S.table.E_orifice_sum + S.table.E_struct_sum; %#ok<AGROW>
                end
            for jj = 1:numel(candCols)
                nm = candCols{jj};
                if ismember(nm, S.table.Properties.VariableNames)
                    val = S.table.(nm);
                        if ~all(isfinite(val(:)))
                            val(~isfinite(val)) = 0;
                        end
                    meta.(nm) = val;
                end
            end
        end

    memo(key) = meta;
    if nargout > 2
        details = S;
    end

end

function T = prepend_baseline_row(T, params, scaled, Opost, lambda, pwr, W)
% PREPEND_BASELINE_ROW Pareto tablosunun başına referans satırı ekler.
% Amaç ve Kapsam: GA ile bulunan çözümleri, başlangıç parametrelerinin
%   yeniden değerlendirilmiş performansı ile kıyaslanabilir kılar.
% Girdiler:
%   T [-]: Pareto ön cephesinden gelen tablo; değişken adları sabit varsayılır.
%   params [-]: Başlangıç parametre yapısı (SI birimleri).
%   scaled [-]: Ölçekli yer hareketi kayıtları seti.
%   Opost [-]: GA sonrası değerlendirme ayarları (thr, penalty).
%   lambda [-]: Penalty ölçek katsayısı.
%   pwr [-]: Penalty güç üssü.
%   W [-]: Penalty ağırlıkları (Δp, Qcap, cav, T, μ bileşenleri).
% Çıktılar:
%   T [-]: İlk satırı baseline değerlendirmesi olan genişletilmiş tablo.
% Varsayımlar ve Sınırlar:
%   - T değişkenleri GA değerlendirmesi ile tutarlıdır.
%   - penalty parametreleri GA sırasında kullanılan değerlerle aynıdır.
%   - run_batch_windowed pencereleme politikaları Opost içinden okunur.
% Ölçüler/Birimler: dP95 [Pa], Qcap95 [-]:, cav_pct [%], T_end [°C], μ_end [Pa·s],
%   E_orifice_sum [J], E_struct_sum [J], energy_tot_sum [J], P_mech_sum [J].
% Yöntem Özeti: params önce encode_params_to_x ile karar vektörüne çevrilir,
%   quant_clamp_x ızgarasına kilitlenir, ardından eval_design_fast ile
%   değerlendirilir. summarize_metrics_table metri̇kleri hesaplayarak baseline
%   satırının hücrelerini doldurur.
    X0 = encode_params_to_x(params);
    X0 = quant_clamp_x(X0);
    [f0, ~, S0] = eval_design_fast(X0, scaled, params, Opost);
    tbl0 = table();
    if isstruct(S0) && isfield(S0,'table') && istable(S0.table)
        tbl0 = S0.table;
    end
    metrics = summarize_metrics_table(tbl0, Opost, lambda, pwr, W);

    T0 = T(1,:);
    T0{1,:} = nan(1,width(T0));

    names = {'d_o_mm','n_orf','Cd0','CdInf','p_exp','Lori_mm','hA_W_perK', ...
             'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref'};
    values = num2cell(X0);
    for idx_loc = 1:numel(names)
        T0 = set_field_safe(T0, names{idx_loc}, values{idx_loc});
    end

    T0 = set_field_safe(T0,'f_pen',f0(1));
    T0 = set_field_safe(T0,'f1',   f0(2));
    T0 = set_field_safe(T0,'f2',   f0(3));
    metric_names = {'PFA_mean','IDR_mean','pen','pen_dP','pen_Qcap','pen_cav', ...
                    'pen_T','pen_mu','x10_max_damperli','a10abs_max_damperli', ...
                    'dP95','Qcap95','cav_pct','T_end','mu_end','Q_q50','Q_q95', ...
                    'dP50','E_orifice_sum','E_struct_sum','energy_tot_sum', ...
                    'E_ratio','P_mech_sum'};
    metric_vals = {metrics.PFA_mean,metrics.IDR_mean,metrics.pen,metrics.pen_dP, ...
                   metrics.pen_Qcap,metrics.pen_cav,metrics.pen_T,metrics.pen_mu, ...
                   metrics.x10_max_damperli,metrics.a10abs_max_damperli,metrics.dP95, ...
                   metrics.Qcap95,metrics.cav_pct,metrics.T_end,metrics.mu_end, ...
                   metrics.Q_q50,metrics.Q_q95,metrics.dP50,metrics.E_orifice_sum, ...
                   metrics.E_struct_sum,metrics.energy_tot_sum,metrics.E_ratio, ...
                   metrics.P_mech_sum};
    for idx_loc = 1:numel(metric_names)
        T0 = set_field_safe(T0, metric_names{idx_loc}, metric_vals{idx_loc});
    end

    T = [T0; T];
end

%% Çıktıların Yazılması ve Raporlama
% Sonuçların dosyalanması ve raporlanması bu bölümde toplanır.
function write_pareto_results(T, outdir)
% WRITE_PARETO_RESULTS Pareto tablosunu kalıcı çıktıya yazar.
% Amaç ve Kapsam: GA sonuçlarını CSV formatında dışa aktararak sürüm takibini
%   ve post-processing adımlarını kolaylaştırır.
% Girdiler:
%   T [-]: Pareto ön cephesine ait tablo (array2table çıktısı).
%   outdir [-]: Yazılacak dizin yolu.
% Çıktılar: Yok; ga_front.csv dosyasını oluşturur.
% Varsayımlar ve Sınırlar:
%   - outdir dizini önceden mevcut olmalıdır.
%   - Tablo sütun adları raporlama araçlarıyla uyumludur.
% Ölçüler/Birimler: Tablo içinde dP95 [Pa], Qcap95 [-]:, cav_pct [%],
%   P_mech_sum [J], energy_tot_sum [J], T_end [°C], μ_end [Pa·s] yer alır.
% Yöntem Özeti: MATLAB writetable çağrısı ile UTF-8 uyumlu CSV üretir.
    writetable(T, fullfile(outdir,'ga_front.csv'));
end


%% Metriklerin Özetlenmesi ve QC
% Enerji, güç ve QC bayrakları bu çekirdek üzerinden tek noktadan yönetilir.
function metrics = summarize_metrics_table(tbl, Opost, lambda, pwr, W)
% SUMMARIZE_METRICS_TABLE pencere metriklerini ceza bileşenleriyle toplar.
% Amaç ve Kapsam: run_batch_windowed çıktısındaki tabloyu okuyarak tasarım
%   bazında performans metriklerini ve pen_raw = Σ W_i·δ_i terimini üretir.
% Girdiler:
%   tbl [-]: run_batch_windowed tarafından üretilen kayıt-tablosu.
%   Opost [-]: QC eşikleri (thr) ve penalty opsiyonlarını içeren yapı.
%   lambda [-]: Penalty ölçek katsayısı.
%   pwr [-]: Penalty güç üssü.
%   W [-]: Penalty ağırlıkları (Δp, Qcap, cav, T, μ bileşenleri).
% Çıktılar:
%   metrics [-]: PFA_mean [m/s²], IDR_mean [-]:, pen [-]:, enerji toplamları [J],
%     P_mech_sum [J], Q kuantilleri [m³/s] ve QC bayraklarını içeren yapı.
% Varsayımlar ve Sınırlar:
%   - tbl sütunları yoksa enforce_finite sıfır doldurması yapılır.
%   - QC eşikleri thr alanında bulunur; eksikse compute_thr ile türetilmiş olmalıdır.
%   - Penalty hesapları boyutsuzdur ve yalnızca pencere bazlı ortalamalara dayanır.
% Ölçüler/Birimler: PFA [m/s²], IDR [-]:, dP95 [Pa], Qcap95 [-]:, cav_pct [%],
%   T_end [°C], μ_end [Pa·s], energy_tot_sum [J], E_orifice_sum [J], E_struct_sum [J],
%   P_mech_sum [J] (Δp·q·Δt integrali ≙ [W·s]).
% Yöntem Özeti: get_numeric_column eksik sütunları yakalar; enforce_finite NaN
%   ve Inf değerlerini temizler. QC bayrakları ok_T, ok_mu, ok_dP, ok_Qcap,
%   ok_cav thr eşikleriyle karşılaştırılır: T_end ≤ T_end_max, μ_end ≥ μ_end_min,
%   dP95 ≤ dP95_max, Qcap95 < Qcap95_max, cav_pct = 0 şartları aranır. Bu bayraklar
%   tasarımın QC başarısını belirler; cav_pct < cav_free olduğunda penalty cezası
%   üretmez ancak ok_cav yalnızca tam sıfır kavitasyonda true olur. Penalty
%   bileşenleri δ_i göreli sapma olarak hesaplanıp pen = λ·(Σ W_i·δ_i)^pwr
%   formunda değerlendirilir.
    if nargin < 1 || isempty(tbl)
        tbl = table();
    end
    if nargin < 2 || isempty(Opost)
        Opost = struct();
    end
    penopt = util_getfield_default(Opost,'penalty', struct());
    cav_free = util_getfield_default(penopt,'cav_free',0.002);
    if nargin < 3 || isempty(lambda)
        lambda = util_getfield_default(penopt,'lambda',10);
    end
    if nargin < 4 || isempty(pwr)
        pwr = util_getfield_default(penopt,'power',2.0);
    end
    if nargin < 5 || isempty(W)
        baseW = util_getfield_default(penopt,'W', struct());
        if ~isstruct(baseW), baseW = struct(); end
        W = struct('dP', util_getfield_default(baseW,'dP',2), ...
                   'Qcap', util_getfield_default(baseW,'Qcap',3), ...
                   'cav', util_getfield_default(baseW,'cav',1), ...
                   'T', util_getfield_default(baseW,'T',0.5), ...
                   'mu', util_getfield_default(baseW,'mu',0.5));
    end
    if ~isstruct(W)
        W = struct();
    end

    rel = @(v,lim) max(0,(v - lim)./max(lim,eps)).^pwr;
    rev = @(v,lim) max(0,(lim - v)./max(lim,eps)).^pwr;

    v_PFA        = enforce_finite(get_numeric_column(tbl,'PFA',0));
    v_IDR        = enforce_finite(get_numeric_column(tbl,'IDR',0));
    v_x10        = enforce_finite(get_numeric_column(tbl,'x10_max_damperli',0));
    v_a10        = enforce_finite(get_numeric_column(tbl,'a10abs_max_damperli',0));
    v_dP95       = enforce_finite(get_numeric_column(tbl,'dP95',0));
    v_Qcap       = enforce_finite(get_numeric_column(tbl,'Qcap95',0));
    v_cav        = enforce_finite(get_numeric_column(tbl,'cav_pct',0));
    v_T_end      = enforce_finite(get_numeric_column(tbl,'T_end',0));
    v_mu         = enforce_finite(get_numeric_column(tbl,'mu_end',0));
    v_Q_q50      = enforce_finite(get_numeric_column(tbl,'Q_q50',0));
    v_Q_q95      = enforce_finite(get_numeric_column(tbl,'Q_q95',0));
    v_dP50       = enforce_finite(get_numeric_column(tbl,'dP50',0));
    v_E_orifice  = enforce_finite(get_numeric_column(tbl,'E_orifice_sum',0));
    v_E_struct   = enforce_finite(get_numeric_column(tbl,'E_struct_sum',0));
    v_P_mech     = enforce_finite(get_numeric_column(tbl,'P_mech_sum',0));

    % Enerji [J] ve mekanik güç integrali [W·s≈J] büyüklükleri disipatif/iş
    % girdisi olduğu için nümerik dalgalanmalarla negatif çıkmaları fiziksel
    % değildir; taban düzeltmesi ile ≥0 garantilenir.
    v_E_orifice = max(v_E_orifice, 0);
    v_E_struct  = max(v_E_struct,  0);
    v_P_mech    = max(v_P_mech,    0);

    metrics = struct();
    metrics.PFA_mean            = mean(v_PFA(:));
    metrics.IDR_mean            = mean(v_IDR(:));
    metrics.x10_max_damperli    = max(v_x10(:));
    metrics.a10abs_max_damperli = max(v_a10(:));
    metrics.dP95                = max(v_dP95(:));
    metrics.Qcap95              = max(v_Qcap(:));
    metrics.cav_pct             = max(v_cav(:));
    metrics.T_end               = max(v_T_end(:));
    metrics.mu_end              = min(v_mu(:));
    metrics.Q_q50               = max(v_Q_q50(:));
    metrics.Q_q95               = max(v_Q_q95(:));
    metrics.dP50                = max(v_dP50(:));
    metrics.E_orifice_sum       = sum(v_E_orifice(:));
    metrics.E_struct_sum        = sum(v_E_struct(:));
    metrics.energy_tot_sum      = metrics.E_orifice_sum + metrics.E_struct_sum;
    metrics.E_ratio             = 0;
    if metrics.E_struct_sum > 0
        metrics.E_ratio = metrics.E_orifice_sum / max(metrics.E_struct_sum, eps);
    end
    metrics.P_mech_sum          = sum(v_P_mech(:));

    thr = util_getfield_default(Opost, 'thr', struct());
    W.dP   = util_getfield_default(W,'dP',2);
    W.Qcap = util_getfield_default(W,'Qcap',3);
    W.cav  = util_getfield_default(W,'cav',1);
    W.T    = util_getfield_default(W,'T',0.5);
    W.mu   = util_getfield_default(W,'mu',0.5);
    metrics.pen_dP   = mean(rel(v_dP95(:), util_getfield_default(thr,'dP95_max',inf)));
    metrics.pen_Qcap = mean(rel(v_Qcap(:), util_getfield_default(thr,'Qcap95_max',inf)));
    cav_lim = util_getfield_default(thr,'cav_pct_max',0);
    cav_eff = max(0, v_cav(:) - cav_free);
    if cav_lim <= 0
        metrics.pen_cav = mean(cav_eff.^pwr);
    else
        cav_lim_eff = max(cav_lim - cav_free, eps);
        metrics.pen_cav = mean((cav_eff ./ max(cav_lim_eff, eps)).^pwr);
    end
    metrics.pen_T    = mean(rel(v_T_end(:), util_getfield_default(thr,'T_end_max',inf)));
    metrics.pen_mu   = mean(rev(v_mu(:), util_getfield_default(thr,'mu_end_min',0)));
    pen_raw          = W.dP*metrics.pen_dP + W.Qcap*metrics.pen_Qcap + ...
                       W.cav*metrics.pen_cav + W.T*metrics.pen_T + W.mu*metrics.pen_mu;
    metrics.pen_raw  = pen_raw;
    metrics.pen      = lambda * pen_raw;
end

function sig = dataset_signature(scaled)
% DATASET_SIGNATURE pencere değerlendirme cache'i için imza üretir.
% Amaç ve Kapsam: scaled kayıt seti değiştiğinde memoize edilmiş sonuçları
%   ayırt etmek; GA boyunca tutarlı veri kullanımını garanti etmek.
% Girdiler:
%   scaled [-]: Ölçekli kayıt yapısı dizisi (name, scale, SaT1 alanları beklenir).
% Çıktılar:
%   sig [-]: jsonencode ile uyumlu string imza.
% Varsayımlar ve Sınırlar:
%   - scaled elemanlarında 'name' ve 'scale' alanları bulunur.
%   - İmza yalnızca temel alanları içerir; ayrıntılı dalga formları dahil edilmez.
% Ölçüler/Birimler: scale [-]:, SaT1 [m/s²].
% Yöntem Özeti: structfun ile alanlar çıkarılır, jsonencode ile deterministik
%   imza oluşturulur.
    if nargin < 1 || isempty(scaled)
        sig = 0;
        return;
    end
    IM = arrayfun(@(s) util_getfield_default(s,'IM',0), scaled);
    PGA = arrayfun(@(s) util_getfield_default(s,'PGA',0), scaled);
    sig = sum(double(IM(:))) + sum(double(PGA(:))) + numel(scaled);
end

function vec = enforce_finite(vec)
% ENFORCE_FINITE NaN veya Inf değerleri sıfırla değiştirir.
% Amaç ve Kapsam: Penalty hesaplamalarında kararlı davranış sağlamak.
% Girdiler:
%   vec [-]: Sayısal vektör veya dizi.
% Çıktılar:
%   vec [-]: NaN/Inf yerine 0 atanmış dizi.
% Varsayımlar ve Sınırlar:
%   - Giriş realdir; kompleks değerler beklenmez.
% Ölçüler/Birimler: Giriş hangi birimde ise çıktı da aynı birimdedir.
% Yöntem Özeti: isfinite maskesi ile eleme yapılır, seçilen elemanlar 0 yapılır.
    if isempty(vec)
        vec = 0;
    end
    if ~isnumeric(vec)
        vec = 0;
        return;
    end
    vec(~isfinite(vec)) = 0;
end

function arr = get_numeric_column(tbl, varName, defaultVal)
% GET_NUMERIC_COLUMN Tablo sütunlarını güvenli şekilde çıkarır.
% Amaç ve Kapsam: Eksik sütunları varsayılan değerle doldurarak ceza ve QC
%   hesaplarında tutarlılık sağlar.
% Girdiler:
%   tbl [-]: MATLAB tablosu.
%   varName [-]: Sütun adı (char veya string).
%   defaultVal [-]: Varsayılan skaler değer.
% Çıktılar:
%   arr [-]: Sayısal sütun vektörü.
% Varsayımlar ve Sınırlar:
%   - defaultVal sayısal olmalıdır.
%   - Sütun mevcut değilse sabit değer döndürülür.
% Ölçüler/Birimler: Sütun birimleri korunur; defaultVal aynı birimde olmalıdır.
% Yöntem Özeti: ismember ile sütun varlığı kontrol edilir, yoksa defaultVal
%   ile dolu vektör üretilir.
    if nargin < 3
        defaultVal = 0;
    end
    if istable(tbl) && ismember(varName, tbl.Properties.VariableNames)
        arr = tbl.(varName);
        if isempty(arr)
            arr = defaultVal;
        end
    else
        arr = defaultVal;
    end
    if ~isnumeric(arr)
        arr = defaultVal;
    end
end

function Trow = set_field_safe(Trow, name, val)
% SET_FIELD_SAFE Tablo satırında güvenli alan ataması yapar.
% Amaç ve Kapsam: Tablo sütunu mevcutsa değer atamak, değilse sessizce geçmek.
% Girdiler:
%   Trow [-]: Tek satırlı tablo.
%   name [-]: Sütun adı.
%   val [-]: Atanacak değer.
% Çıktılar:
%   Trow [-]: Güncellenmiş tablo satırı.
% Varsayımlar ve Sınırlar:
%   - Trow tek satırlı olmalıdır.
% Ölçüler/Birimler: Giriş değeri hangi birimde ise aynı kalır.
% Yöntem Özeti: ismember ile sütun varlığını kontrol eder ve {} ataması kullanır.
    if ~istable(Trow)
        return;
    end
    if ~ismember(name, Trow.Properties.VariableNames)
        return;
    end
    if iscell(Trow.(name))
        Trow.(name) = {val};
    else
        if isnumeric(val)
            Trow.(name) = val;
        else
            Trow.(name) = NaN;
        end
    end
end

function x = encode_params_to_x(params)
% ENCODE_PARAMS_TO_X Parametre yapısını karar vektörüne dönüştürür.
% Amaç ve Kapsam: GA karar uzayındaki sıralamayı koruyarak params alanlarını
%   [d_o_mm, n_orf, Cd0, CdInf, p_exp, Lori_mm, hA_W_perK, Dp_mm, d_w_mm,
%   D_m_mm, n_turn, mu_ref] vektörüne dönüştürmek.
% Girdiler:
%   params [-]: build_params tarafından oluşturulan yapı.
% Çıktılar:
%   x [-]: Karar vektörü; geometrik terimler [mm], katsayılar boyutsuz, mu_ref [-]:.
% Varsayımlar ve Sınırlar:
%   - params içinde gerekli alanlar mevcuttur.
% Ölçüler/Birimler: mm değerleri orifis çapları ve uzunluklarına aittir.
% Yöntem Özeti: util_getfield_default ile eksik alanlara NaN atanır ve vektör
%   sıraya göre doldurulur.
    orf = util_getfield_default(params,'orf', struct());
    thermal = util_getfield_default(params,'thermal', struct());

    x = nan(1,12);
    x(1) = 1e3 * util_getfield_default(orf,'d_o',3.0e-3);
    x(2) = util_getfield_default(params,'n_orf',6);
    x(3) = util_getfield_default(orf,'Cd0',0.61);
    x(4) = util_getfield_default(orf,'CdInf',0.80);
    x(5) = util_getfield_default(orf,'p_exp',1.10);
    x(6) = 1e3 * util_getfield_default(params,'Lori',0.15);
    x(7) = util_getfield_default(thermal,'hA_W_perK',600);
    x(8) = 1e3 * util_getfield_default(params,'Dp',0.150);
    x(9) = 1e3 * util_getfield_default(params,'d_w',0.012);
    x(10) = 1e3 * util_getfield_default(params,'D_m',0.120);
    x(11) = util_getfield_default(params,'n_turn',8);
    x(12) = util_getfield_default(params,'mu_ref',0.9);

    if ~isfinite(x(2))
        x(2) = 6;
    end
end




%% Tasarım Değişkenlerinin Izgaraya Kilitlenmesi
% Karar vektörü bileşenleri belirlenmiş adımlara ve sınırlara sıkıştırılır.
function xq = quant_clamp_x(x)
% QUANT_CLAMP_X tasarım vektörünü ızgara adımlarına uyarlar.
% Amaç ve Kapsam: GA popülasyonunda yapılandırılmış aramayı korumak için her
%   bileşeni önceden tanımlı adımlarla kuantize eder.
% Girdiler:
%   x [-]: [d_o_mm, n_orf, Cd0, CdInf, p_exp, Lori_mm, hA_W_perK, Dp_mm,
%        d_w_mm, D_m_mm, n_turn, mu_ref] karar vektörü.
% Çıktılar:
%   xq [-]: Aynı boyutta kuantize edilmiş vektör.
% Varsayımlar ve Sınırlar:
%   - x bileşenleri alt/üst sınırlar arasında olmalıdır.
%   - Tamsayı bileşenler (n_orf, n_turn) en yakın tam değere yuvarlanır.
% Ölçüler/Birimler: mm değerleri mm olarak kalır; katsayılar boyutsuz.
% Yöntem Özeti: util_quantize_step adımları uygular ve clamp ile sınırlar.
%   Izgara adımları tablo halinde: d_o_mm 0.05 [mm], Cd 0.01 [-]:, p_exp 0.05 [-]:,
%   Lori_mm 1 [mm], hA_W_perK 25 [W/K], Dp_mm 1 [mm], d_w_mm 0.5 [mm],
%   D_m_mm 5 [mm], mu_ref 0.05 [Pa·s (ref)]; n_orf ve n_turn tam sayı olarak
%   yuvarlanır.
    % quant_clamp_x — Karar vektörünü GA/eval ile aynı ızgaraya kilitler
    % Izgara adımları: d_o_mm 0.05 [mm], Cd0/CdInf 0.01 [-]:, p_exp 0.05 [-]:
    %   Lori_mm 1 [mm], hA_W_perK 25 [W/K], Dp_mm 1 [mm], d_w_mm 0.5 [mm],
    %   D_m_mm 5 [mm], mu_ref 0.05 [Pa·s referans]. n_orf, n_turn tam sayı ≥1.
    % Boyut analizi: [mm] → karar uzayı, dönüştürme decode_params_from_x içinde m'ye yapılır.
    xq = x;
    if isvector(x), xq = x(:)'; end
    xq(:,1) = util_quantize_step(xq(:,1),0.05);
    xq(:,3) = util_quantize_step(xq(:,3),0.01);
    xq(:,4) = util_quantize_step(xq(:,4),0.01);
    xq(:,5) = util_quantize_step(xq(:,5),0.05);
    if size(xq,2) >= 6,  xq(:,6)  = util_quantize_step(xq(:,6),1);   end
    if size(xq,2) >= 7,  xq(:,7)  = util_quantize_step(xq(:,7),25);  end
    if size(xq,2) >= 8,  xq(:,8)  = util_quantize_step(xq(:,8),1);   end
    if size(xq,2) >= 9,  xq(:,9)  = util_quantize_step(xq(:,9),0.5); end
    if size(xq,2) >=10,  xq(:,10) = util_quantize_step(xq(:,10),5);  end
    if size(xq,2) >=11,  xq(:,11) = round(xq(:,11));                end
    if size(xq,2) >=12,  xq(:,12) = util_quantize_step(xq(:,12),0.05); end
    xq(:,2) = round(max(xq(:,2),1));
    if size(xq,2) >=11, xq(:,11) = round(max(xq(:,11),1)); end
    if isvector(x), xq = xq(:)'; end
end


%% Parametre Kod Çözümü (mm → m) ve Model Kurulumu
% mm ile verilen ızgara değişkenlerini SI tabanlı modele geri yazar.
function P = decode_params_from_x(params_base_, x_)
% DECODE_PARAMS_FROM_X karar vektörünü parametre yapısına çözer.
% Amaç ve Kapsam: GA çıktısını build_params tarafından beklenen SI birimli
%   yapı formatına döndürmek.
% Girdiler:
%   params_base_ [-]: Referans parametre yapısı.
%   x_ [-]: [d_o_mm, n_orf, Cd0, CdInf, p_exp, Lori_mm, hA_W_perK, Dp_mm,
%          d_w_mm, D_m_mm, n_turn, mu_ref] karar vektörü.
% Çıktılar:
%   P [-]: SI birimlerine dönüştürülmüş parametre yapısı.
% Varsayımlar ve Sınırlar:
%   - Geometri girişleri mm cinsindedir ve 1e-3 ile çarpılarak metreye çevrilir.
%   - n_orf ve n_turn tamsayıya yuvarlanır.
% Ölçüler/Birimler: d_o = d_o_mm·1e−3 [m], Lori = Lori_mm·1e−3 [m], Dp_mm·1e−3 [m],
%   d_w_mm·1e−3 [m], D_m_mm·1e−3 [m], hA_W_perK [W/K], mu_ref [Pa·s].
% Yöntem Özeti: params_base_ kopyalanır, 1e−3 faktörleri ile mm→m dönüşümü
%   yapılır ve build_params çağrısı ile bağlı parametreler güncellenir.
    % decode_params_from_x — Karar vektöründen fiziksel parametreleri türetir
    % Birim dönüşümleri: [mm] → [m], Cd [-]:, p_exp [-]:, hA_W_perK [W/K], μ_ref [Pa·s]
    % Boyut kontrolü: Uzunluklar 1e-3 çarpanı ile metreye çevrilir.
    d_o_mm = x_(1); n_orf = round(x_(2));
    Cd0 = x_(3); CdInf = x_(4); p_exp = x_(5);
    Lori_mm = x_(6); hA_W_perK = x_(7);
    Dp_mm = x_(8); d_w_mm = x_(9); D_m_mm = x_(10);
    n_turn = round(x_(11)); mu_ref = x_(12);
    P = params_base_;
    P.orf.d_o = d_o_mm * 1e-3;         % mm'den m'ye
    P.n_orf   = n_orf;
    P.orf.Cd0   = Cd0;
    P.orf.CdInf = CdInf;
    P.orf.p_exp = p_exp;
    P.Lori   = Lori_mm * 1e-3;
    P.mu_ref = mu_ref;
    P.Dp     = Dp_mm * 1e-3;
    P.d_w    = d_w_mm * 1e-3;
    P.D_m    = D_m_mm * 1e-3;
    P.n_turn = n_turn;
    if isfield(P,'thermal')
        P.thermal.hA_W_perK = hA_W_perK;
    end
    P = build_params(P);
end
function [state, options, optchanged] = ga_out_best_pen(options, state, flag, scaled, params, optsEval)
% GA_OUT_BEST_PEN GA iterasyonlarında en iyi tasarımı raporlar.
% Amaç ve Kapsam: Iteratif pen, f1, f2 değerlerini stdout'a yazmak ve
%   meta veriyi güncel tutmak.
% Girdiler:
%   options [-]: GA seçenek yapısı.
%   state [-]: GA durum yapısı (Population, Score vb.).
%   flag [-]: GA callback bayrağı.
%   scaled [-]: Ölçekli kayıt seti (memo cache için).
%   params [-]: Parametre yapısı.
%   optsEval [-]: Değerlendirme seçenekleri.
% Çıktılar:
%   state [-]: Değiştirilmeden geri döner.
%   options [-]: Değiştirilmeden geri döner.
%   optchanged [-]: false; callback ayarları değiştirmez.
% Varsayımlar ve Sınırlar:
%   - state.Score sütunu [f_pen, f1, f2] sıralamasını içerir.
% Ölçüler/Birimler: f_pen [-]:, f1 [m/s²], f2 [-]:.
% Yöntem Özeti: Score matrisini leksikografik olarak sıralar, en iyi bireyin
%   eval_design_fast değerini hesaplar ve fprintf ile kaydeder.
    optchanged = false;
    if strcmp(flag,'iter') && ~isempty(state.Score)
        scores = state.Score;                 % [f_pen, f1, f2]
        [~, order] = sortrows(scores, [1 2 3]);
        idx = order(1);
        bestx = state.Population(idx,:);
        [f_curr, ~] = eval_design_fast(bestx, scaled, params, optsEval);
        fprintf('Gen %d: pen=%g f1=%g f2=%g\\n', state.Generation, f_curr(1), f_curr(2), f_curr(3));
    end
end

function [scaled, params, T1] = prepare_inputs(optsGA)
% PREPARE_INPUTS GA sürücüsü için veri ve parametreleri toparlar.
% Amaç ve Kapsam: Varsayılan parametreleri kurar ve yer hareketi kayıtlarını
%   yükler.
% Girdiler:
%   optsGA [-]: Yükleme seçenekleri (load_opts alanı isteğe bağlı).
% Çıktılar:
%   scaled [-]: Ölçekli kayıt seti.
%   params [-]: Parametre yapısı (SI birimleri).
%   T1 [s]: Birinci mod periyodu.
% Varsayımlar ve Sınırlar:
%   - load_opts yoksa parametreler_defaults tarafından belirlenen kayıtlar kullanılır.
% Ölçüler/Birimler: T1 [s], diğer alanlar build_params birimlerine uygundur.
% Yöntem Özeti: parametreler_defaults çağrılır, load_opts varsa load_ground_motions
%   ile veri seti okunur.

if nargin < 1 || isempty(optsGA)
    optsGA = struct();
end

[params, T1] = parametreler_defaults();
load_opts = util_getfield_default(optsGA, 'load_opts', []);
if isempty(load_opts)
    [~, scaled] = load_ground_motions(T1);
else
    [~, scaled] = load_ground_motions(T1, load_opts);
end

end

function [params, T1] = parametreler_defaults()
% PARAMETRELER_DEFAULTS yapı ve damper varsayılanlarını kurar.
% Amaç ve Kapsam: 10 katlı kesme çerçevesi ve damper özelliklerini SI
%   birimlerinde oluşturmak.
% Girdiler: Yok.
% Çıktılar:
%   params [-]: Yapı ve damper parametrelerini içeren yapı.
%   T1 [s]: Birinci mod periyodu.
% Varsayımlar ve Sınırlar:
%   - Kat sayısı n = 10 sabittir.
%   - Kütle, rijitlik ve sönüm değerleri eş dağıtılmıştır.
% Ölçüler/Birimler: m [kg], k [N/m], c [N·s/m], uzunluklar [m], sıcaklık [°C],
%   viskozite [Pa·s], yoğunluk [kg/m³].
% Yöntem Özeti: Kat kütleleri, rijitlik ve sönüm matrisleri oluşturulur,
%   eig ile T1 hesaplanır, damper geometri verileri mm→m olarak tanımlanır,
%   orifis ve termal parametreler eklenir.

%% --- 1) Structure (10 story shear building) ---
n  = 10;
m  = 360e3 * ones(n,1);
k  = 6.5e8 * ones(n,1);
c  = 6.2e6 * ones(n,1);
story_height = 3.0;

M  = diag(m);
K  = zeros(n);
C0 = zeros(n);
for i = 1:n
    kL = k(i);    cL = c(i);
    if i < n
        kU = k(i+1); cU = c(i+1);
    else
        kU = 0;    cU = 0;
    end
    K(i,i)   = kL + kU;   C0(i,i)   = cL + cU;
    if i > 1
        K(i,i-1) = -kL;  C0(i,i-1) = -cL;
    end
    if i < n
        K(i,i+1) = -kU;  C0(i,i+1) = -cU;
    end
end

[~,D] = eig(K,M);
w  = sqrt(sort(diag(D),'ascend'));
T1 = 2*pi/w(1);

%% --- 2) Damper geometry and material data ---
Dp     = 0.125;
Lgap   = 0.055;
d_o    = 3.0e-3;
Lori   = 0.10;
mu_ref = 0.9;

Kd     = 1.6e9;
Ebody  = 2.1e11;
Gsh    = 79e9;
d_w    = 12e-3;
D_m    = 80e-3;
n_turn = 8;

%% --- 3) Orifice and thermal parameters ---
rho   = 850;
n_orf = 6;

orf = struct('Cd0', 0.61, ...
             'CdInf', 0.8, ...
             'Rec', 3000, ...
             'p_exp', 1, ...
             'p_amb', 1e5, ...
             'p_cav_eff', -1e7, ...
             'cav_sf', 1, ...
             'd_o', d_o, ...
             'veps', 0.1);

T0_C    = 25;
T_ref_C = 25;
b_mu    = -0.013;
thermal = struct();
thermal.hA_W_perK = 450;
thermal.T_env_C   = 25;
thermal.max_iter  = 3;
thermal.tol_K     = 0.5;
thermal.relax     = 0.5;
thermal.dT_max    = 80;

steel_to_oil_mass_ratio = 1.5;
n_dampers_per_story    = 1;
story_mask             = ones(n-1,1);
cp_oil   = 1800;
cp_steel = 500;
resFactor = 12;

c_lam_cap      = 2e7;
c_lam_min_frac = 0.6;
c_lam_min_abs  = 1e5;

    cfg = struct();
    cfg.on = struct('pressure_force', false, 'mu_floor', false, 'thermal_feedback', true);
cfg.compat_simple = false;
cfg.num = struct('softmin_eps', 1e6, 'mu_min_phys', 0.6, 'dP_cap', 2e7, 'Qcap_scale', 1);

params = struct('M',M,'C0',C0,'K',K,'Dp',Dp,'Lgap',Lgap,'d_w',d_w, ...
    'D_m',D_m,'n_turn',n_turn,'mu_ref',mu_ref,'Lori',Lori,'Kd',Kd, ...
    'Ebody',Ebody,'Gsh',Gsh,'rho',rho,'n_orf',n_orf,'orf',orf, ...
    'thermal',thermal,'T0_C',T0_C,'T_ref_C',T_ref_C,'b_mu',b_mu, ...
    'c_lam_cap',c_lam_cap,'c_lam_min_frac',c_lam_min_frac, ...
    'c_lam_min_abs',c_lam_min_abs,'cp_oil',cp_oil,'cp_steel',cp_steel, ...
    'steel_to_oil_mass_ratio',steel_to_oil_mass_ratio, ...
    'n_dampers_per_story',n_dampers_per_story,'story_mask',story_mask, ...
    'resFactor',resFactor,'cfg',cfg,'story_height',story_height);

params = build_params(params);

end

%% Kayıtların Pencereli Analizi (Arias)
% Arias yoğunluğu pencereleri (t5–t95) üzerinden toplu değerlendirme yapılır.
function [summary, all_out] = run_batch_windowed(scaled, params, opts)
% RUN_BATCH_WINDOWED pencere tabanlı çoklu kayıt çözümleyicisidir.
% Amaç ve Kapsam: scaled kayıtlarının her biri için run_one_record_windowed
%   çağırarak pencere metriklerini ve QC bayraklarını toparlar.
% Girdiler:
%   scaled [-]: Ölçekli kayıt dizisi.
%   params [-]: Yapı ve damper parametreleri (build_params çıktısı).
%   opts [-]: Pencere, termal politika ve QC eşiklerini içeren yapı.
% Çıktılar:
%   summary [-]: Tablo ve bayraklardan oluşan özet yapı.
%   all_out [-]: Her kayıt için ayrıntılı çıktı hücre dizisi.
% Varsayımlar ve Sınırlar:
%   - params.thermal.hA_W_perK mevcut olmalıdır.
%   - opts.thr eksikse compute_thr ile türetilir.
% Ölçüler/Birimler: Zaman [s], PFA [m/s²], IDR [-]:, T_end [°C], μ_end [Pa·s],
%   P_mech_sum [J], Q [m³/s].
% Yöntem Özeti: rbw_prepare_inputs dizileri hazırlar; rbw_record_loop Arias
%   yoğunluğu pencerelerini kullanarak run_one_record_windowed çağrılarını
%   sıralar. Elde edilen vars yapısı rbw_build_summary_table ile QC ve ceza
%   özetine dönüştürülür.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'thr'), opts.thr = struct(); end
opts.thr = compute_thr(params, opts.thr);

% İstenirse kayıtların yürütülme sırası yeniden düzenlenir
if isfield(opts,'order_perm')
    scaled = scaled(opts.order_perm);
end

assert(isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), ...
    'run_batch_windowed: params.thermal.hA_W_perK eksik');

%% Girdi Hazırlığı
n = numel(scaled);
vars = rbw_prepare_inputs(n, params, opts);

%% Kayıt Döngüsü
vars = rbw_record_loop(scaled, params, opts, vars);

%% Özet Tablo
summary = rbw_build_summary_table(vars, opts);
all_out = vars.all_out;

end

function vars = rbw_prepare_inputs(n, params, opts)
% RBW_PREPARE_INPUTS pencere analizi için veri kapsayıcılarını oluşturur.
% Amaç ve Kapsam: run_batch_windowed içindeki kayıt döngüsüne giriş sağlayacak
%   dizileri başlatır.
% Girdiler:
%   n [-]: Kayıt sayısı.
%   params [-]: Parametre yapısı (termal politika bilgisi için).
%   opts [-]: Termal ve sıralama politikası seçenekleri.
% Çıktılar:
%   vars [-]: Boş diziler ve meta bilgileri içeren yapı.
% Varsayımlar ve Sınırlar:
%   - opts.thermal_reset ∈ {'each','carry','cooldown'} olarak kullanılır.
% Ölçüler/Birimler: Zaman [s], sıcaklık [°C], viskozite [Pa·s], basınç [Pa].
% Yöntem Özeti: util_getfield_default ile politikalar okunur, her metrik için
%   sıfır vektörler oluşturulur.
vars = struct();
vars.all_out = cell(n,1);
vars.names    = cell(n,1);
vars.scale    = zeros(n,1);
vars.SaT1     = zeros(n,1);
vars.t5       = zeros(n,1);
vars.t95      = zeros(n,1);
vars.coverage = zeros(n,1);

policy_val = util_getfield_default(opts,'thermal_reset','each');
order_val  = util_getfield_default(opts,'order','natural');
vars.policy_col = repmat({policy_val}, n,1);
vars.order_col  = repmat({order_val}, n,1);
if isfield(opts,'cooldown_s')
    cooldown_val = opts.cooldown_s;
else
    cooldown_val = NaN;
end
vars.cooldown_col = repmat(cooldown_val, n,1);

vars.PFA    = zeros(n,1);
vars.IDR    = zeros(n,1);
vars.dP95   = zeros(n,1);
vars.Qcap95 = zeros(n,1);
vars.cav_pct = zeros(n,1);
vars.zeta1_hot       = zeros(n,1);
vars.z2_over_z1_hot  = zeros(n,1);
vars.P_mech_sum      = zeros(n,1);
vars.Re_max      = zeros(n,1);
vars.Q_q95  = zeros(n,1);
vars.Q_q50  = zeros(n,1);
vars.dP50   = zeros(n,1);
vars.x10_max_damperli = zeros(n,1);
vars.a10abs_max_damperli = zeros(n,1);
vars.E_orifice_sum = zeros(n,1);
vars.E_struct_sum  = zeros(n,1);
vars.energy_tot_sum = zeros(n,1);
vars.E_ratio   = zeros(n,1);
vars.qc_pass   = false(n,1);

vars.T_start    = zeros(n,1);
vars.T_end      = zeros(n,1);
vars.mu_end     = zeros(n,1);
vars.clamp_hits = zeros(n,1);

vars.Dp_mm_col   = repmat(util_getfield_default(params,'Dp_mm',NaN), n,1);
vars.mu_ref_col  = repmat(util_getfield_default(params,'mu_ref',NaN), n,1);

end

function vars = rbw_record_loop(scaled, params, opts, vars)
% RBW_RECORD_LOOP her kaydı Arias penceresi ile işler.
% Amaç ve Kapsam: run_one_record_windowed çağrıları yaparak vars yapılarını
%   doldurur.
% Girdiler:
%   scaled [-]: Ölçekli kayıt dizisi.
%   params [-]: Parametre yapısı.
%   opts [-]: Termal ve pencere seçenekleri.
%   vars [-]: rbw_prepare_inputs tarafından başlatılmış yapı.
% Çıktılar:
%   vars [-]: Metrikler ve meta bilgilerle güncellenmiş yapı.
% Varsayımlar ve Sınırlar:
%   - run_one_record_windowed termal politikayı prev_ts üzerinden takip eder.
% Ölçüler/Birimler: Zaman [s], PFA [m/s²], IDR [-]:, dP95 [Pa], Qcap95 [-]:,
%   cav_pct [%], T_end [°C], μ_end [Pa·s], P_mech_sum [J].
% Yöntem Özeti: Önceki kaydın termal durumunu taşımak için prev_ts değişkeni
%   saklanır; çıkış metrikleri vars alanlarına aktarılır.
prev_ts = [];
for k = 1:numel(scaled)
    rec = scaled(k);
    out = run_one_record_windowed(rec, params, opts, prev_ts);
    prev_ts = out.ts;
    vars.all_out{k} = out; %#ok<AGROW>

    vars.names{k}    = out.name;
    vars.scale(k)    = out.scale;
    vars.SaT1(k)     = out.SaT1;
    vars.t5(k)       = out.win.t5;
    vars.t95(k)      = out.win.t95;
    vars.coverage(k) = out.win.coverage;

    vars.T_start(k)    = out.T_start;
    vars.T_end(k)      = out.T_end;
    vars.mu_end(k)     = out.mu_end;
    vars.clamp_hits(k) = out.clamp_hits;
    m_nom = out.metr;
    vars.PFA(k)    = m_nom.PFA;
    vars.IDR(k)    = m_nom.IDR;
    vars.dP95(k)   = m_nom.dP95;
    vars.Qcap95(k) = m_nom.Qcap95;
    vars.cav_pct(k)= m_nom.cav_pct;
    vars.zeta1_hot(k)       = util_getfield_default(m_nom,'zeta1_hot',NaN);
    vars.z2_over_z1_hot(k)  = util_getfield_default(m_nom,'z2_over_z1_hot',NaN);
    vars.P_mech_sum(k)      = util_getfield_default(m_nom,'P_mech_sum',NaN);
    vars.Re_max(k)          = util_getfield_default(m_nom,'Re_max',NaN);
    vars.Q_q95(k)  = util_getfield_default(m_nom,'Q_q95',NaN);
    vars.Q_q50(k)  = util_getfield_default(m_nom,'Q_q50',NaN);
    vars.dP50(k)   = util_getfield_default(m_nom,'dP50',NaN);
    vars.x10_max_damperli(k) = util_getfield_default(m_nom,'x10_max_damperli',NaN);
    vars.a10abs_max_damperli(k) = util_getfield_default(m_nom,'a10abs_max_damperli',NaN);
    vars.E_orifice_sum(k) = util_getfield_default(m_nom,'E_orifice_sum',NaN);
    vars.E_struct_sum(k)  = util_getfield_default(m_nom,'E_struct_sum',NaN);
    vars.energy_tot_sum(k) = util_getfield_default(m_nom,'energy_tot_sum',NaN);
    vars.E_ratio(k)   = util_getfield_default(m_nom,'E_ratio',NaN);
    vars.qc_pass(k)   = out.qc_pass;

end
end

function summary = rbw_build_summary_table(vars, opts)
% RBW_BUILD_SUMMARY_TABLE kayıt bazlı metriklerden özet yapı oluşturur.
% Amaç ve Kapsam: vars alanlarını tabloya çevirerek QC bayraklarını hesaplar.
% Girdiler:
%   vars [-]: rbw_record_loop tarafından doldurulan yapı.
%   opts [-]: QC eşiklerini içeren yapı (thr alanı gereklidir).
% Çıktılar:
%   summary [-]: table ve all_out alanlarını içeren özet.
% Varsayımlar ve Sınırlar:
%   - opts.thr alanı T_end_max [°C], μ_end_min [Pa·s], dP95_max [Pa],
%     Qcap95_max [-]: değerlerini sağlar.
% Ölçüler/Birimler: Tablo sütunları PFA [m/s²], IDR [-]:, dP95 [Pa],
%   Qcap95 [-]:, cav_pct [%], enerji [J], T_end [°C], μ_end [Pa·s].
% Yöntem Özeti: QC bayrakları thr ile karşılaştırılır; ok_cav yalnızca cav_pct = 0
%   olduğunda true olup sıfır kavitasyon şartını temsil eder. Ancak penalty
%   tarafında cav_free ölü-bölgesi cav_pct < cav_free için ceza üretmez.
%   qc_reason alanı başarısız eşikleri listeler, summary.table tabloya dönüştürür.
summary = struct();

thr = opts.thr;
ok_T    = vars.T_end   <= thr.T_end_max;
ok_mu   = vars.mu_end  >= thr.mu_end_min;
ok_dP   = vars.dP95    <= thr.dP95_max;
ok_Qcap = vars.Qcap95  <  thr.Qcap95_max;
ok_cav  = vars.cav_pct == 0;
qc_reason = strings(numel(vars.names),1);
for r = 1:numel(vars.names)
    bad = {};
    if ~ok_T(r),    bad{end+1}='T';  end %#ok<AGROW>
    if ~ok_mu(r),   bad{end+1}='mu'; end %#ok<AGROW>
    if ~ok_dP(r),   bad{end+1}='dP'; end %#ok<AGROW>
    if ~ok_Qcap(r), bad{end+1}='Qcap'; end %#ok<AGROW>
    if ~ok_cav(r),  bad{end+1}='cav'; end %#ok<AGROW>
    qc_reason(r) = strjoin(bad,',');
end

summary.table = table(vars.names, vars.scale, vars.SaT1, vars.t5, vars.t95, vars.coverage, vars.policy_col, vars.order_col, vars.cooldown_col, ...
    vars.PFA, vars.IDR, vars.dP95, vars.Qcap95, vars.cav_pct, vars.zeta1_hot, vars.z2_over_z1_hot, vars.P_mech_sum, vars.Re_max, ...
    vars.Q_q95, vars.Q_q50, vars.dP50, vars.x10_max_damperli, vars.a10abs_max_damperli, vars.E_orifice_sum, vars.E_struct_sum, vars.energy_tot_sum, vars.E_ratio, vars.qc_pass, ...
    vars.T_start, vars.T_end, vars.mu_end, vars.clamp_hits, vars.Dp_mm_col, vars.mu_ref_col, ...
    ok_T, ok_mu, ok_dP, ok_Qcap, ok_cav, qc_reason, ...
    'VariableNames', {'name','scale','SaT1','t5','t95','coverage','policy','order','cooldown_s', ...
    'PFA','IDR','dP95','Qcap95','cav_pct','zeta1_hot','z2_over_z1_hot','P_mech_sum','Re_max', ...
    'Q_q95','Q_q50','dP50','x10_max_damperli','a10abs_max_damperli','E_orifice_sum','E_struct_sum','energy_tot_sum','E_ratio','qc_pass', ...
    'T_start','T_end','mu_end','clamp_hits','Dp_mm','mu_ref','ok_T','ok_mu','ok_dP','ok_Qcap','ok_cav','qc_reason'});

summary.all_out = vars.all_out;
end

%% Kayıt Analiz Fonksiyonu
function out = run_one_record_windowed(rec, params, opts, prev_ts)
% RUN_ONE_RECORD_WINDOWED tek bir kayıt için Arias pencereli analizi yürütür.
% Amaç ve Kapsam: mck_with_damper çözümü ile zaman serilerini üretir ve
%   compute_metrics_windowed üzerinden pencere metriklerini toplar.
% Girdiler:
%   rec [-]: Ölçekli yer hareketi kaydı (alanlar: name, scale, SaT1, ag, dt).
%   params [-]: Yapı ve damper parametreleri.
%   opts [-]: Pencere ayarları (window), termal politika (thermal_reset, cooldown_s)
%          ve QC seçenekleri.
%   prev_ts [-]: Önceki kayıt sonundan taşınan termal zaman serisi (opsiyonel).
% Çıktılar:
%   out [-]: name, scale, win, metr, ts, qc_pass, T_start [°C], T_end [°C],
%         μ_end [Pa·s], clamp_hits [-], enerji ve akış metriklerini içeren yapı.
% Varsayımlar ve Sınırlar:
%   - rec.ag zaman serisi [m/s²] olarak sağlanır.
%   - Termal reset politikası 'each', 'carry' veya 'cooldown' değerlerinden biridir.
% Ölçüler/Birimler: Zaman [s], hız [m/s], ivme [m/s²], sıcaklık [°C],
%   viskozite [Pa·s], basınç [Pa], akış [m³/s], enerji [J].
% Yöntem Özeti: util_make_arias_window ile t5–t95 aralığı oluşturulur;
%   mck_with_damper kayıt boyunca çözülür. P_mech ≈ Σ(Δp·q·Δt) ilişkisi üzerinden
%   compute_metrics_windowed enerji ve güç metriklerini hesaplar. Termal modlar
%   cooldown_s süresi kadar bekleme, doğrudan taşıma veya her kayıt sonrası reset
%   seçenekleri ile uygulanır.

% varsayılan argümanlar
if nargin < 4, prev_ts = []; end
if nargin < 3 || isempty(opts), opts = struct(); end

if isfield(opts,'thermal_reset') && strcmpi(opts.thermal_reset,'cooldown')
    if ~isfield(opts,'cooldown_s') || isempty(opts.cooldown_s) || isnan(opts.cooldown_s)
        opts.cooldown_s = 60;
    end
    opts.cooldown_s = max(opts.cooldown_s,0);
end

% QC eşikleri (eksik alanlar varsayılanlarla doldurulur)
if ~isfield(opts,'thr'), opts.thr = struct(); end
thr = compute_thr(params, opts.thr);

assert(isfield(params,'thermal') && isfield(params.thermal,'hA_W_perK'), ...
    'run_one_record_windowed: params.thermal.hA_W_perK eksik');


%% Pencere Hazırlığı
if isfield(opts,'window') && ~isempty(opts.window)
    wfields = fieldnames(opts.window);
    wargs = cell(1,2*numel(wfields));
    for k = 1:numel(wfields)
        wargs{2*k-1} = wfields{k};
        wargs{2*k}   = opts.window.(wfields{k});
    end
    win = util_make_arias_window(rec.t, rec.ag, wargs{:});
else
    win = util_make_arias_window(rec.t, rec.ag);
end


% ham ve ölçekli kayıt karşılaştırması kaldırıldı (üst seviyede kullanılmıyor)

%% Termal Sıfırlama
Tinit = params.T0_C;
if isfield(opts,'thermal_reset')
    mode = opts.thermal_reset;
else
    mode = 'each';
end

switch mode
    case 'each'
        Tinit = params.T0_C;
    case 'carry'
        if ~isempty(prev_ts) && isfield(prev_ts,'T_oil')
            Tinit = prev_ts.T_oil(end);
        else
            Tinit = params.T0_C;
        end
    case 'cooldown'
        if ~isempty(prev_ts) && isfield(prev_ts,'T_oil')
            Tprev = prev_ts.T_oil(end);
        else
            Tprev = params.T0_C;
        end
        if isfield(opts,'cooldown_s')
            td = opts.cooldown_s;
        else
            td = 0;
        end
        C_th = util_compute_Cth_effective(params);
        hA = params.thermal.hA_W_perK;
        Tenv = params.thermal.T_env_C;
        Tinit = Tenv + (Tprev - Tenv) * exp(-hA*td / C_th);
    otherwise
        Tinit = params.T0_C;
end

%% Damperli Çözüm
% Damperli çözüm mck_with_damper fonksiyonu üzerinden yürütülür.

mu_ref_eff   = params.mu_ref;
c_lam0_eff   = params.c_lam0;
if isfield(params,'cfg') && isstruct(params.cfg) && isfield(params.cfg,'on') && isstruct(params.cfg.on) && ...
        isfield(params.cfg.on,'mu_floor') && params.cfg.on.mu_floor
    mu_min_phys = NaN;
    if isfield(params.cfg,'num') && isstruct(params.cfg.num) && ...
            isfield(params.cfg.num,'mu_min_phys') && isfinite(params.cfg.num.mu_min_phys)
        mu_min_phys = params.cfg.num.mu_min_phys;
    end
    if isfinite(mu_min_phys) && (mu_ref_eff < mu_min_phys)
        mu_ref_eff = mu_min_phys;
        scale_mu = mu_ref_eff / max(params.mu_ref, eps);
        c_lam0_eff = params.c_lam0 * scale_mu;
    end
end
    [x,a_rel,ts] = mck_with_damper(rec.t, rec.ag, params.M, params.C0, params.K, ...
    params.k_sd, c_lam0_eff, params.Lori, params.orf, params.rho, params.Ap, ...
    params.Ao, params.Qcap_big, mu_ref_eff, params.thermal, Tinit, ...
    params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, params.Lgap, ...
    params.cp_oil, params.cp_steel, params.steel_to_oil_mass_ratio, ...
    params.story_mask, params.n_dampers_per_story, params.resFactor, params.cfg);

metr = compute_metrics_windowed(rec.t, x, a_rel, rec.ag, ts, params.story_height, win, params);

% Nominal koşu için son termal ve viskozite durumları kaydedilir
T_start = Tinit;
if isfield(ts,'T_oil')
    T_end = ts.T_oil(end);
    Tmax = params.T0_C + params.thermal.dT_max;
    clamp_hits = sum(diff(ts.T_oil >= Tmax) > 0);
else
    T_end = NaN;
    clamp_hits = NaN;
end
if isfield(ts,'mu')
    mu_end = ts.mu(end);
else
    mu_end = NaN;
end

qc_pass = (metr.cav_pct <= thr.cav_pct_max) && ...
          (metr.dP95 <= thr.dP95_max) && ...
          (metr.Qcap95 <= thr.Qcap95_max) && ...
          (T_end <= thr.T_end_max) && ...
          (mu_end >= thr.mu_end_min);

%% Çıktıların Derlenmesi
out = struct('name', rec.name, ...
    'scale', util_getfield_default(rec,'scale',1), ...
    'SaT1', util_getfield_default(rec,'IM',NaN), ...
    'win', win, 'metr', metr, 'ts', ts, 'qc_pass', qc_pass, ...
    'T_start', T_start, 'T_end', T_end, 'mu_end', mu_end, ...
    'clamp_hits', clamp_hits, ...
    'PFA', metr.PFA, 'IDR', metr.IDR, 'dP95', metr.dP95, 'Qcap95', metr.Qcap95, ...
    'cav_pct', metr.cav_pct, 't5', win.t5, 't95', win.t95, 'coverage', win.coverage);

% Yeni parametreleri isteğe bağlı olarak kaydet
param_fields = {'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref'};
for ii = 1:numel(param_fields)
    fn = param_fields{ii};
    if isfield(params, fn)
        out.(fn) = params.(fn);
    end
end
end

function metr = compute_metrics_windowed(t, x, a_rel, ag, ts, story_height, win, params)
% COMPUTE_METRICS_WINDOWED pencere içi performans metriklerini hesaplar.
% Amaç ve Kapsam: Arias penceresiyle sınırlanan veri üzerinde hidrolik ve yapısal
%   performans göstergelerini üretmek.
% Girdiler:
%   t [s]: Zaman vektörü.
%   x [m]: Kat göreli yer değiştirmeleri matrisi.
%   a_rel [m/s²]: Göreli ivmeler.
%   ag [m/s²]: Taban ivmesi.
%   ts [-]: Termal ve hidrolik zaman serileri (dP_orf [Pa], Q [m³/s], mu [Pa·s], dvel [m/s], cav_mask [-]:).
%   story_height [m]: Kat yüksekliği.
%   win [-]: Arias penceresi (idx mantıksal maskesi, t5, t95, coverage).
%   params [-]: Model parametreleri (Qcap_big [m³/s], Ap [m²], Ao [m²], rho [kg/m³], orf bilgileri).
% Çıktılar:
%   metr [-]: PFA [m/s²], IDR [-]:, dP95 [Pa], Qcap95 [-]:, cav_pct [%], enerji
%           bileşenleri [J], P_mech_sum [J], Re_max [-]: gibi alanları içeren yapı.
% Varsayımlar ve Sınırlar:
%   - win.idx pencere içinde true/false maskesidir.
%   - ts alanları pencerede yeterli uzunluğa sahiptir.
% Ölçüler/Birimler: Enerji terimleri [J], güç terimleri [W], kavite oranı [%],
%   Reynolds sayısı [-]:, sıcaklık [°C] (ts içinde), viskozite [Pa·s].
% Yöntem Özeti: Tepe ivme ve ötelenme hesaplanır; quantile ile dP ve Q
%   istatistikleri türetilir. P_mech ≈ Σ(Δp·q·Δt) ilişkisi kullanılarak enerji
%   bileşenleri elde edilir; yapısal enerjiye göre E_ratio hesaplanır. Kavitasyon
%   maskesi ile cav_pct belirlenir.

idx = win.idx;

%% Temel Tepki
% Bu bölüm, tepe kat ivmesi ve katlar arası göreli hareketler gibi yapısal
% yanıtın temel göstergelerini hesaplar.
% Tepe katın mutlak ivmesi
a_top_abs = a_rel(:,end) + ag(:);
metr.PFA = max(abs(a_top_abs(idx)));

% Pencere içindeki tepe kat mutlak yerdeğiştirmesi ve ivmesi
metr.x10_max_damperli    = max(abs(x(idx,end)));
metr.a10abs_max_damperli = max(abs(a_top_abs(idx)));

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
metr.Re_max = NaN;
if all(isfield(params,{'Qcap_big','Ap','Ao','rho','orf'})) && isfield(ts,'mu') && isfield(ts,'dvel') && ...
        isstruct(params.orf) && all(isfield(params.orf,{'veps','d_o'}))
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
end

% 50. ve 95. yüzdelik değerler
dP_q50         = quantile(abs_dP, 0.50);
dP_q95         = quantile(abs_dP, 0.95);
Q_q50          = quantile(abs_Q, 0.50);
Q_q95          = quantile(abs_Q, 0.95);
Qcap95 = quantile(Qcap_ratio, 0.95);
story_force_q95= quantile(abs_story_force, 0.95);

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
metr.E_orifice_sum = ts.E_orf(end);
metr.E_struct_sum  = ts.E_struct(end);
metr.E_ratio   = metr.E_orifice_sum / max(metr.E_struct_sum, eps);
% Toplam enerji ve ortalama mekanik güç
metr.energy_tot_sum = metr.E_orifice_sum + metr.E_struct_sum;
metr.P_mech_sum = NaN;
if isfield(ts,'P_sum') && ~isempty(ts.P_sum)
    metr.P_mech_sum = mean(ts.P_sum(idx));
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


% Parametrelerden gerekli alanlar kayıt altına alınır
param_fields = {'Dp_mm','d_w_mm','D_m_mm','n_turn','mu_ref'};
for ii = 1:numel(param_fields)
    fn = param_fields{ii};
    if isfield(params, fn)
        metr.(fn) = params.(fn);
    end
end

end

% MCK_WITH_DAMPER damperli çok serbestlik dereceli sistemi çözer.
% Amaç ve Kapsam: Yapısal M, C, K matrislerini viskoz damper modeliyle
%   birleştirir, termal ve hidrolik durumları günceller ve zaman serilerini üretir.
% Girdiler:
%   t [s]: Zaman vektörü.
%   ag [m/s²]: Taban ivmesi.
%   M [kg]: Kütle matrisi.
%   C [N·s/m]: Sönüm matrisi (pasif + laminer katkı).
%   K [N/m]: Rijitlik matrisi.
%   k_sd [N/m]: Seri yay rijitliği.
%   c_lam0 [N·s/m]: Laminer sönüm referansı.
%   Lori [m]: Çalışma uzunluğu.
%   orf [-]: Orifis parametre yapısı (Cd0, CdInf, p_exp, vs.).
%   rho [kg/m³], Ap [m²], Ao [m²], Qcap [m³/s], mu_ref [Pa·s]: Akış ve viskozite parametreleri.
%   thermal [-]: Termal ayarlar (hA_W_perK [W/K], T_env_C [°C], vb.).
%   T0_C [°C], T_ref_C [°C], b_mu [1/°C], c_lam_min [N·s/m], c_lam_cap [N·s/m], Lgap [m].
%   cp_oil [J/(kg·K)], cp_steel [J/(kg·K)], steel_to_oil_mass_ratio [-]:, story_mask [-]:.
%   n_dampers_per_story [-]: Damper adedi; resFactor [-]: çözünürlük ölçeği; cfg [-]: çözüm ayarları.
% Çıktılar:
%   x [m]: Kat göreli yer değiştirmeleri.
%   a_rel [m/s²]: Kat göreli ivmeleri.
%   ts [-]: Basınç, debi, enerji ve termal tarihçeleri içeren yapı.
% Varsayımlar ve Sınırlar:
%   - Laminer sönüm c_lam ≥ c_lam_min ve ≤ c_lam_cap ile sınırlandırılır.
%   - Orifis modeli Cd(Re) bağıntısını ve kavitasyon sınırını uygular.
% Ölçüler/Birimler: Δp [Pa], q [m³/s], enerji [J], sıcaklık [°C], viskozite [Pa·s].
% Yöntem Özeti: Doğrusal MCK sistemi laminer ve seri yay katkılarıyla güncellenir,
%   orifis basınç düşümü ve kavitasyon `util_softmin` ile harmanlanır, termal
%   durum `util_compute_Cth_effective` ile izlenir ve zaman integrasyonu `ode15s`
%   kullanılarak gerçekleştirilir.
function [x,a_rel,ts] = mck_with_damper(t,ag,M,C,K, k_sd,c_lam0,Lori, orf,rho,Ap,Ao,Qcap, mu_ref, ...
    thermal, T0_C,T_ref_C,b_mu, c_lam_min,c_lam_cap,Lgap, ...
    cp_oil,cp_steel, steel_to_oil_mass_ratio, story_mask, ...
    n_dampers_per_story, resFactor, cfg)
% mck_with_damper — Doğrusal MCK + doğrusal-yay (k_sd) + laminer (c_lam) + orifis (Δp_orf) damper modeli
% Yöneten denklem:
%   M * ẍ + C0 * ẋ + K * x + f_damper = - M * r * a_g
% Kat i için damper kuvveti:
%   f_i = k_sd * Δx_i  +  c_lam * Δv_i  +  A_p * Δp_orf,i * sgn(Δv_i)
%   [k_sd]: N/m  |  [c_lam]: N·s/m  |  [Δx]: m  |  [Δv]: m/s  |  [A_p]: m²  |  [Δp_orf]: Pa
% Orifis akışı ve basınç düşümü:
%   q = Q_cap * tanh( (A_p / Q_cap) * sqrt(Δv^2 + ε_v^2) )             [q]: m³/s
%   Re = (ρ * q * d_o) / (A_o * μ)                                      [ρ]: kg/m³, [μ]: Pa·s, [d_o]: m, [A_o]: m²
%   C_d(Re) = C_d,∞ - (C_d,∞ - C_d,0) / (1 + (Re / Re_c)^p_exp)
%   Δp_kv = ½ ρ ( q / (C_d A_o) )²
%   Δp_cav = max( (p_up - p_cav,eff) * cav_sf, 0 ),   p_up = p_amb + |F_lin| / A_p
%   Δp_orf = softmin(Δp_kv, Δp_cav; ε)   (p-state yok, kararlı harmanlama)
% Enerji denklik kontrolü (pencere içi):
%   Σ_i f_i Δv_i  =  Σ_i (c_lam Δv_i²)  +  Σ_i (Δp_orf,i q_i)  +  d/dt(E_kin + E_pot)
% Boyut kontrolü:
%   [Δp q] = (Pa)(m³/s) = N·m/s = W,  [c_lam Δv²] = (N·s/m)(m²/s²) = N·m/s = W
% Çoklu damper ölçeği:
%   multi = (#damper/story) · (aktiflik maskesi); bu faktör laminer ve orifis kuvvetlerine de uygulanır.
% Termo-viskozite geri beslemesi şu an tanısal; dinamiğe bağlanmadı.
%% Girdi Parametreleri
    n = size(M,1); r = ones(n,1);
    thermal_feedback_on = true;
    if isfield(cfg,'on') && isfield(cfg.on,'thermal_feedback')
        thermal_feedback_on = logical(cfg.on.thermal_feedback);
    end %#ok<NASGU>
    % Termal geri besleme kancası; ileride viskozite-dinamik bağlaması için kullanılacak.
    agf = griddedInterpolant(t,ag,'linear','nearest');
    z0 = zeros(2*n,1);
    opts= odeset('RelTol',1e-3,'AbsTol',1e-6);

    % Kat vektörleri
    nStories = n-1;
    mask = story_mask(:);  if numel(mask)==1,  mask  = mask *ones(nStories,1); end
    ndps = n_dampers_per_story(:); if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
    multi = (mask .* ndps).';
    multi_row = multi(:).';
    Nvec = 1:nStories; Mvec = 2:n;

    % Başlangıç sıcaklığı ve viskozitesi
    Tser = T0_C*ones(numel(t),1);
    mu_abs = mu_ref;
    c_lam = c_lam0;

%% ODE Çözümü
    % Faz 6: dev_force i��in Qcap ve softmin parametrelerini �nceden ayarla
    Qcap_eff = Qcap;
    if isfield(cfg,'num') && isfield(cfg.num,'Qcap_scale') && isfinite(cfg.num.Qcap_scale)
        Qcap_eff = max(1e-9, Qcap * cfg.num.Qcap_scale);
    end
    orf_loc = orf;
    if isfield(cfg,'num') && isfield(cfg.num,'softmin_eps') && isfinite(cfg.num.softmin_eps)
        orf_loc.softmin_eps = cfg.num.softmin_eps;
    end

    odef = @(tt,z) [ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - dev_force(tt,z(1:n),z(n+1:end),c_lam,mu_abs) - M*r*agf(tt) ) ];
    sol  = ode15s(odef,[t(1) t(end)],z0,opts);
    z    = deval(sol,t).';
    x    = z(:,1:n); v = z(:,n+1:end);

    drift = x(:,Mvec) - x(:,Nvec);
    dvel  = v(:,Mvec) - v(:,Nvec);
    % Faz 3: Lineer parca sadece yay (laminer viskoz katkı tarafında)
    F_lin = k_sd*drift;

    % Faz 6: Qcap ölçeği ve softmin eps opsiyonu
    % NOTE: Qcap_eff and orf_loc are pre-initialized above for dev_force.
    params = struct('Ap',Ap,'Qcap',Qcap_eff,'orf',orf_loc,'rho',rho,...
                    'Ao',Ao,'mu',mu_abs,'F_lin',F_lin,'Lori',Lori);
    [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force(dvel, params);
    % Ek diagnostikler: dP_kv ve dP_cav (kv ve kavitasyon limitleri)
    qmag_loc = Qcap_eff * tanh( (Ap/Qcap_eff) * sqrt(dvel.^2 + orf.veps^2) );
    Re_loc   = (rho .* qmag_loc .* max(orf.d_o,1e-9)) ./ max(Ao*mu_abs,1e-9);
    Cd_loc0  = orf.Cd0; Cd_locInf = orf.CdInf; Rec_loc = orf.Rec; pexp_loc = orf.p_exp;
    Cd_loc = Cd_locInf - (Cd_locInf - Cd_loc0) ./ (1 + (Re_loc./max(Rec_loc,1)).^pexp_loc);
    Cd_loc = max(min(Cd_loc,1.2),0.2);
    dP_kv_loc = 0.5*rho .* ( qmag_loc ./ max(Cd_loc.*Ao,1e-12) ).^2;
    p_up_loc  = orf.p_amb + abs(F_lin)./max(Ap,1e-12);
    dP_cav_loc= max( (p_up_loc - orf.p_cav_eff).*orf.cav_sf, 0 );
    % --- Kavitasyon maskesi (fiziksel eşiklere dayalı) ---
    % Tolerans: ölçüm/hesap dalgalanmalarını yutmak için %5 dP_cav veya min 1e4 Pa
    tol_loc = max(1e4, 0.05*max(dP_cav_loc,[],'all'));  % [Pa]
    Qcap_ratio_loc = abs(Q) ./ max(Qcap_eff, eps);       % [-]
    cav_mask_loc = (dP_kv_loc > (dP_cav_loc - tol_loc)) | (Qcap_ratio_loc > 0.98);
    % Paketleme: cav_mask_loc [-] mantıksal; Δp_kv ≳ Δp_cav veya Q/Qcap → 1 bölgeleri işaretlenir.
    % ts.cav_mask = cav_mask_loc;   % boyut: [Nt x nStories], mantıksal
    F_p = F_lin + F_orf;

    multi_mat = repmat(multi_row, size(dvel,1), 1);
    F_visc_per = c_lam .* dvel;
    F_visc_total = F_visc_per .* multi_mat;
    F_orf_total = F_orf .* multi_mat;


    % Geometri ölçeklendirmesi R sadece montajda uygulanır
    F_story = F_lin + F_visc_total + F_orf_total;
    % NOT: Enerji tarafında P_visc_per = c_lam .* (dvel.^2) [W] ve P_orf_per = Δp_orf .* |Q| [W],
    % ve toplamlar multi_mat ile ölçeklendiğinden, dinamik (F·Δv) ile enerji (P_toplam) aynı ölçekle tutarlı hale geldi.
    P_visc_per = c_lam .* (dvel.^2);
    P_sum = sum( (P_visc_per + P_orf_per) .* multi_mat, 2 );
    P_orf_tot = sum(P_orf_per .* multi_mat, 2);
    % Yapısal güç kat toplam kuvvetini kullanır; ekstra çarpan kullanılmaz
    P_struct_tot = sum(F_story .* dvel, 2);
    E_orf = cumtrapz(t, P_orf_tot);
    E_struct = cumtrapz(t, P_struct_tot);

%% Termal Hesap (Phase 7: iki-düğüm diagnostik; dinamiğe geri besleme yok)
    nDtot = sum(multi);
    V_oil_per = resFactor*(Ap*(2*Lgap));
    m_oil_tot = nDtot*(rho*V_oil_per);
    m_steel_tot = steel_to_oil_mass_ratio*m_oil_tot;
    C_oil   = max(m_oil_tot*cp_oil,   eps);
    C_steel = max(m_steel_tot*cp_steel, eps);
    T_o = Tser; T_s = T0_C*ones(numel(t),1);
    hA_os   = util_getfield_default(thermal, 'hA_os',    thermal.hA_W_perK);
    hA_o_env= util_getfield_default(thermal, 'hA_o_env', thermal.hA_W_perK);
    hA_s_env= util_getfield_default(thermal, 'hA_s_env', thermal.hA_W_perK);
    dtv = diff(t);
    for k=1:numel(t)-1
        Pk = 0.5*(P_sum(k)+P_sum(k+1));
        dT_o = ( Pk - hA_os*(T_o(k)-T_s(k)) - hA_o_env*(T_o(k)-thermal.T_env_C) ) / C_oil;
        dT_s = ( + hA_os*(T_o(k)-T_s(k)) - hA_s_env*(T_s(k)-thermal.T_env_C) ) / C_steel;
        T_o(k+1) = T_o(k) + dtv(k)*dT_o;
        T_s(k+1) = T_s(k) + dtv(k)*dT_s;
        T_o(k+1) = min(max(T_o(k+1), T0_C), T0_C + thermal.dT_max);
        T_s(k+1) = min(max(T_s(k+1), T0_C), T0_C + thermal.dT_max);
    end
    mu = mu_ref*exp(b_mu*(T_o - T_ref_C));

%% Çıktı Hesabı
    % İvme için düğüm kuvvetleri
    F = zeros(numel(t),n);
    F(:,Nvec) = F(:,Nvec) - F_story;
    F(:,Mvec) = F(:,Mvec) + F_story;
    a_rel = ( -(M\(C*v.' + K*x.' + F.')).' - ag.*r.' );

    ts = struct('dvel', dvel, 'story_force', F_story, 'Q', Q, ...
        'dP_orf', dP_orf, 'cav_mask', cav_mask_loc, 'P_sum', P_sum, ...
        'E_orf', E_orf, 'E_struct', E_struct, 'T_oil', T_o, 'mu', mu, 'c_lam', c_lam);

%% İç Fonksiyonlar
    function Fd = dev_force(tt,x_,v_,c_lam_loc,mu_abs_loc)
        drift_ = x_(Mvec) - x_(Nvec);
        dvel_  = v_(Mvec) - v_(Nvec);
        multi_row   = multi(:).';              % [1 x nStories]
        multi_col   = multi_row(:);            % [nStories x 1]
        % Sütun yönelimli etkin parametreler
        % Faz 3: Lineer parçada sadece yay
        F_lin_ = k_sd * drift_;
        F_lin_ = F_lin_(:);
        params = struct('Ap',Ap,'Qcap',Qcap_eff,'orf',orf_loc,'rho',rho,...
                        'Ao',Ao,'mu',mu_abs_loc,'F_lin',F_lin_,'Lori',Lori);
        [F_orf_, ~, ~, ~] = calc_orifice_force(dvel_, params);
        F_orf_ = F_orf_(:);
        dvel_  = dvel_(:);
        F_orf_eff_  = F_orf_(:)  .* multi_col;       % orifis kuvveti + multi
        F_visc_raw_ = c_lam_loc .* dvel_;            % laminer kuvvet
        F_visc_raw_ = F_visc_raw_(:);
        F_visc_eff_ = F_visc_raw_ .* multi_col;      % laminer + multi
        F_story_    = F_lin_ + F_visc_eff_ + F_orf_eff_;
        F_story_ = F_story_(:);
        Fd = zeros(n,1);
        Fd(Nvec) = Fd(Nvec) - F_story_;
        Fd(Mvec) = Fd(Mvec) + F_story_;
    end
    function [F_orf, dP_orf, Q, P_orf_per] = calc_orifice_force(dvel, params)
        % calc_orifice_force — Orifis damper akış ve basınç düşümü modeli
        % Tanımlar: q [m³/s], Δp_kv [Pa], Δp_cav [Pa], Δp_orf [Pa], F_orf [N], P_orf [W]
        % P_orifis = Δp_orf · |q| [W]  (laminer güç c_lam · Δv² ayrı tutulur; çifte sayım yok)
        % C_d(Re) korelasyonu ve sınırları (0.2 ≤ C_d ≤ 1.2) ve Re → 0/∞ limitleri
        % Enerji denklik notu: Δp_orf · q boyut kontrolü ile N·m/s = W sağlar.
        % Boyut analizi: Δp [Pa], q [m³/s], F_orf [N], P_orf [W] ⇒ 1 Pa = 1 N/m².
        % Varsayımlar: p-state yok, tanh saturasyonu (ε_v) küçük hızlarda nümerik güvenlik sağlar.

        % Saturated volumetric flow magnitude (stability)
        qmag = params.Qcap * tanh( (params.Ap/params.Qcap) * sqrt(dvel.^2 + params.orf.veps^2) );

        % Reynolds and discharge coefficient (clamped)
        Re   = (params.rho .* qmag .* max(params.orf.d_o,1e-9)) ./ max(params.Ao*params.mu,1e-9);
        Cd0   = params.orf.Cd0;
        CdInf = params.orf.CdInf;
        p_exp = params.orf.p_exp;
        Rec   = params.orf.Rec;
        Cd    = CdInf - (CdInf - Cd0) ./ (1 + (Re./max(Rec,1)).^p_exp);
        Cd    = max(min(Cd, 1.2), 0.2);

        % kv-only drop
        dP_kv  = 0.5*params.rho .* ( qmag ./ max(Cd.*params.Ao,1e-12) ).^2;

        % Cavitation soft-limit via softmin
        p_up   = params.orf.p_amb + abs(params.F_lin)./max(params.Ap,1e-12);
        dP_cav = max( (p_up - params.orf.p_cav_eff).*params.orf.cav_sf, 0 );
        epsm = 1e5;
        if isfield(params,'orf') && isfield(params.orf,'softmin_eps') && isfinite(params.orf.softmin_eps)
            epsm = params.orf.softmin_eps;
        end
        dP_orf = util_softmin(dP_kv, dP_cav, epsm);

        % Force sign from velocity (no p-states)
        sgn = dvel ./ sqrt(dvel.^2 + params.orf.veps^2);
        F_orf = dP_orf .* params.Ap .* sgn;

        % Diagnostics (positive)
        Q = qmag;
        P_orf_per = dP_orf .* abs(Q);
    end
end

% LOAD_GROUND_MOTIONS fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   T1 [-]: SI uyumlu giriş parametresi.
%   opts [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   records [-]: Hesaplanmış çıktı.
%   scaled [-]: Hesaplanmış çıktı.
%   meta [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function [records, scaled, meta] = load_ground_motions(T1, opts)
%LOAD_GROUND_MOTIONS Zemin hareketi kayitlarini yukler, on isler ve olcekler.
%   [RAW, SCALED, META] = LOAD_GROUND_MOTIONS(T1, OPTS) fonksiyonu
%   acc_matrix.mat dosyasindaki acc_matrix1, acc_matrix2, ... degiskenlerini
%   okuyarak her kayit icin zaman vektoru, ivme ve temel buyuklukleri
%   iceren bir yapi dizisi dondurur.
%
%   T1 girilirse kayitlarin 5%% sonumlu T1 periyodundaki veya tanimli
%   bir periyot bandindaki yapay spektral ivmeleri (PSA) hesaplanir ve
%   tum kayitlar ortak bir hedef IM degerine olceklendirilir.
%
%   OPTS alanlari:
%       hp_cut   - yuksek gecis filtresi frekansi [Hz]
%       IM_mode  - 'T1' veya 'band'
%       band_fac - T1 cevresindeki band carpani [alt ust]
%       band_N   - bant icindeki periyot sayisi
%       s_bounds - olcek katsayisi sinirlari [min max]
%       verbose  - islemler hakkinda bilgi mesajlari (true/false)
%
%   Ciktilar:
%       records - ham kayitlar
%       scaled  - olceklenmis kayitlar
%       meta    - islemlere dair ozet bilgiler

%% Girdi Parametreleri
if nargin < 2, opts = struct(); end
if ~isfield(opts,'verbose'), opts.verbose = true; end
hp_cut   = util_getfield_default(opts,'hp_cut',0.05);   % yuksek gecis [Hz]
IM_mode  = util_getfield_default(opts,'IM_mode','band');
band_fac = util_getfield_default(opts,'band_fac',[0.8 1.2]);
band_N   = util_getfield_default(opts,'band_N',21);
s_bounds = util_getfield_default(opts,'s_bounds',[0.2 2.2]);

%% MAT dosyasını Yükle
raw = load('acc_matrix.mat');
fn  = fieldnames(raw);
records = struct('name',{},'t',{},'ag',{},'dt',{},'duration',{},'PGA',{},'PGV',{},'IM',{},'scale',{});

%% Ön İşleme
for k = 1:numel(fn)
    A = raw.(fn{k});
    t  = A(:,1);   ag = A(:,2);

    % 3A) Zaman vektorunu duzenle
    [t,iu] = unique(t,'stable'); ag = ag(iu);
    dt = median(diff(t));
    assert(max(abs(diff(t) - dt)) < 1e-6, 'Zaman örnekleme aralığı düzensiz');
    t  = (t(1):dt:t(end)).';
    ag = interp1(A(:,1), A(:,2), t, 'linear');

    % 3B) Temel birim ve filtreleme kontrolleri
    assert(max(abs(ag)) < 100, 'İvme büyüklüğü birim hatasına işaret ediyor');
    ag = detrend(ag,0);             % ortalama
    ag = detrend(ag,1);             % dogrusal trend
    try
        Fs = 1/dt;
        Wn = hp_cut/(Fs/2);
        [b,a] = butter(2,Wn,'high');
        ag = filtfilt(b,a,ag);
    catch
        % Butter/Filtfilt yoksa yalnizca detrend uygulanir
    end

    % 3C) Metaveri
    duration = t(end) - t(1);
    v  = cumtrapz(t,ag);
    PGA = max(abs(ag));
    PGV = max(abs(v));

    records(end+1) = struct('name',fn{k},'t',t,'ag',ag,'dt',dt, ...
                             'duration',duration,'PGA',PGA,'PGV',PGV, ...
                             'IM',[],'scale',1); %#ok<AGROW>
end

%% Yükleme Özeti
if opts.verbose, fprintf('Toplam %d zemin hareketi kaydı yüklendi\\n', numel(records)); end
for k = 1:numel(records)
    r = records(k);
    if opts.verbose, fprintf('%2d) %-12s dt=%6.4f s dur=%6.2f s PGA=%7.3f PGV=%7.3f\\n', ...
        k, r.name, r.dt, r.duration, r.PGA, r.PGV); end
end

scaled = [];
meta = struct();

%% IM Hesabı ve Ölçekleme (T1 sağlandığında)
if nargin >= 1 && ~isempty(T1)
    %% IM Hesabı
    for k = 1:numel(records)
        records(k).IM = compute_IM(records(k).t, records(k).ag, IM_mode, T1, band_fac, band_N);
    end

    %% Hedef IM Seçimi ve TRIM
    IM = [records.IM];
    IM_low  = max(s_bounds(1)*IM);
    IM_high = min(s_bounds(2)*IM);
    targetIM0 = median(IM);
    targetIM  = min(max(targetIM0, IM_low), IM_high);
    doClip = (IM_low > IM_high);

    max_iter = 3; dropped = {};
    for it = 1:max_iter
        IM_low  = max(s_bounds(1)*IM);
        IM_high = min(s_bounds(2)*IM);
        if IM_low <= IM_high, break; end
        [~,idx] = max(abs(log(IM) - median(log(IM))));
        dropped{end+1} = records(idx).name; %#ok<AGROW>
        records(idx) = [];
        IM(idx) = [];
    end
    if ~isempty(dropped) && opts.verbose, fprintf('TRIM: ayıklanan uç değerler = %s\n', strjoin(dropped,', ')); end
    IM_low  = max(s_bounds(1)*IM);
    IM_high = min(s_bounds(2)*IM);
    targetIM0 = median(IM);
    targetIM  = min(max(targetIM0, IM_low), IM_high);
    doClip = (IM_low > IM_high);

    %% Ölçekleme
    scaled = records; n_clipped = 0; s_all = zeros(1,numel(records));
    for k = 1:numel(records)
        s_raw = targetIM / records(k).IM;
        if doClip
            s = min(max(s_raw, s_bounds(1)), s_bounds(2));
            if abs(s - s_raw) > 1e-12, n_clipped = n_clipped + 1; end
        else
            s = s_raw;
        end
        s_all(k) = s;

        scaled(k).ag    = s * records(k).ag;
        scaled(k).PGA   = s * records(k).PGA;
        scaled(k).PGV   = s * records(k).PGV;
        scaled(k).scale = s;
        scaled(k).s_clipped = doClip && (abs(s - s_raw) > 1e-12);
        scaled(k).trimmed   = false;

        scaled(k).IM = compute_IM(scaled(k).t, scaled(k).ag, IM_mode, T1, band_fac, band_N);
    end

    %% Hata ve Log
    err = abs([scaled.IM] - targetIM) / max(targetIM, eps) * 100;
    if strcmpi(IM_mode,'band')
        modeStr = 'band';
    else
        modeStr = 'PSA@T1';
    end
    clipCount = n_clipped * doClip;
    if opts.verbose, fprintf('Hedef IM = %.3f (%s). Maks hata = %.2f%% | uygun aralık=[%.3f, %.3f] | s_min=%.2f s_max=%.2f | KIRPILAN=%d\\n', ...
        targetIM, modeStr, max(err), IM_low, IM_high, min(s_all), max(s_all), clipCount); end

    meta = struct('IM_mode', IM_mode, 'band_fac', band_fac, 's_bounds', s_bounds);
    if exist('dropped','var'), meta.TRIM_names = dropped; else, meta.TRIM_names = {}; end
end
end

%% ==== Yerel Fonksiyonlar ====
% COMPUTE_IM fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   t [-]: SI uyumlu giriş parametresi.
%   ag [-]: SI uyumlu giriş parametresi.
%   mode [-]: SI uyumlu giriş parametresi.
%   T1 [-]: SI uyumlu giriş parametresi.
%   band_fac [-]: SI uyumlu giriş parametresi.
%   band_N [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   IM [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function IM = compute_IM(t, ag, mode, T1, band_fac, band_N)
%COMPUTE_IM Kayit icin hedef IM degerini hesaplar.
%   IM = COMPUTE_IM(T, AG, MODE, T1, BAND_FAC, BAND_N) fonksiyonu,
%   verilen zaman vektoru t ve ivme ag icin, secilen modda
%   yapay spektral ivme (IM) dondurur.

zeta = 0.05;
if strcmpi(mode,'band')
    Tgrid = linspace(band_fac(1)*T1, band_fac(2)*T1, band_N);
    Sa = zeros(size(Tgrid));
    for i = 1:numel(Tgrid)
        Sa(i) = calc_psa(t, ag, Tgrid(i), zeta);
    end
    IM = exp(mean(log(Sa + eps)));
else
    IM = calc_psa(t, ag, T1, zeta);
end
end

% CALC_PSA fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   t [-]: SI uyumlu giriş parametresi.
%   ag [-]: SI uyumlu giriş parametresi.
%   T [-]: SI uyumlu giriş parametresi.
%   zeta [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   Sa [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function Sa = calc_psa(t, ag, T, zeta)
%CALC_PSA Tek bir kayit icin yapay spektral ivme hesabi.
%   Sa = CALC_PSA(t, ag, T, zeta) fonksiyonu, T periyotlu ve zeta
%   sonumlu tek serbestlik dereceli osilatorun mutlak ivme tepkisinin
%   maksimum degerini dondurur.

w = 2*pi / T;
agf = griddedInterpolant(t, ag, 'linear', 'nearest');

odef = @(tt, y)[ y(2);
                 -2*zeta*w*y(2) - w*w*y(1) - agf(tt) ];

y0 = [0;0];
[~, y] = ode45(odef, t, y0);

x  = y(:,1);
xd = y(:,2);
xdd = -2*zeta*w*xd - w*w*x - ag;
abs_acc = xdd + ag;
Sa = max(abs(abs_acc));
end

% BUILD_PARAMS fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   params [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   params [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function params = build_params(params)
% build_params — Hidrolik ve yapısal türetilmiş alanların hesaplanması
% Ap = π D_p²/4 [m²], A_o = n_orf π d_o²/4 [m²], k_sd [N/m], c_lam0 [N·s/m]
% c_lam0 Poiseuille laminer akışından türetilir: Δp = 128 μ L q /(π d⁴) ⇒
%   F_lam = c_lam0 Δv, c_lam0 ≈ 128 μ L / (π d⁴) · A_p² / L_gap varsayımıyla,
%   paralel kanalların eşdeğer uzunluğu L_gap ve silindirik delik kabulüyle.
% Qcap_big teorik orifis tavanı (Δp_cap [Pa]) ile belirlenir.
% Varsayımlar: Newtonyen yağ (μ sabit referans), sabit sıcaklık T_ref, rijit piston.
% Boyut analizi: Ap,Ao [m²], Qcap_big [m³/s], k_sd [N/m], c_lam0 [N·s/m], Δp_cap [Pa].

if nargin < 1 || ~isstruct(params)
    params = struct();
end

% Reuse existing utility for core damper quantities
params = util_recompute_damper_params(params);

    % Qcap_big — Büyük basınç farkı altında (Δp ≤ dP_cap) teorik azami debi [m³/s]
    % Formül: Qcap_big = C_d,∞ * A_o * sqrt( 2 Δp_cap / ρ )
    % Öneri: Δp_cap ~ 5–30 MPa aralığı (ürün ve güvenlik kısıtlarına göre seçilir). Varsayılan 20 MPa.
    if isfield(params,'orf') && isfield(params.orf,'CdInf') && ...
            isfield(params,'Ao') && isfield(params,'rho')
        dP_cap = util_getfield_default( ...
                     util_getfield_default(params,'cfg',struct()), ...
                     'num', struct());
        dP_cap = util_getfield_default(dP_cap,'dP_cap', 2e7);  % [Pa], varsayılan 20 MPa
        params.Qcap_big = max(params.orf.CdInf * params.Ao, 1e-9) * ...
            sqrt(2*dP_cap / params.rho);
    end

% Minimum laminar damping based on c_lam0
if isfield(params,'c_lam_min_abs') && isfield(params,'c_lam_min_frac') && ...
        isfield(params,'c_lam0')
    params.c_lam_min = max(params.c_lam_min_abs, ...
        params.c_lam_min_frac * params.c_lam0);
end

end

% PARPOOL_HARD_RESET fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   nWorkers [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   status [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function status = parpool_hard_reset(nWorkers)
% PARPOOL_HARD_RESET parpool'u sıfırlayıp iş parçacıklarını kısıtlar.
% Eski işleri temizleyerek güvenli bir havuz açılışı sağlar.

status = true;

if nargin<1 || isempty(nWorkers)
    nWorkers = feature('numcores');
else
    validateattributes(nWorkers, {'numeric'}, {'scalar','integer','positive'});
end

try
    c = parcluster('Processes');

    if ~isempty(c.Jobs)
        delete(c.Jobs);  % çökmüş işleri temizle
    end

    p = gcp('nocreate');
    if isempty(p) || ~isvalid(p)
        parpool(c, min(nWorkers, c.NumWorkers));
    end
catch ME
    warning('parpool_hard_reset:Failed', 'Parallel pool could not be opened: %s', ME.message);
    status = false;
    return;
end

pctRunOnAll maxNumCompThreads(1);          % CPU aşırı kullanımını engelle
pctRunOnAll set(0,'DefaultFigureVisible','off');
pctRunOnAll setenv('OMP_NUM_THREADS','1');
pctRunOnAll setenv('MKL_NUM_THREADS','1');

end

% UTIL_SOFTMIN fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   a [-]: SI uyumlu giriş parametresi.
%   b [-]: SI uyumlu giriş parametresi.
%   epsm [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   y [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function y = util_softmin(a, b, epsm)
%UTIL_SOFTMIN Smooth minimum used for cavitation blending.
y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
end

% UTIL_RECOMPUTE_DAMPER_PARAMS fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   params [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   params [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function params = util_recompute_damper_params(params)
%UTIL_RECOMPUTE_DAMPER_PARAMS Refresh derived damper quantities.
if ~isstruct(params), return; end

if isfield(params,'Dp_mm'),    params.Dp    = params.Dp_mm/1000; end
if isfield(params,'d_w_mm'),   params.d_w   = params.d_w_mm/1000; end
if isfield(params,'D_m_mm'),   params.D_m   = params.D_m_mm/1000; end
if isfield(params,'Lori_mm'),  params.Lori  = params.Lori_mm/1000; end
if isfield(params,'orf') && isfield(params.orf,'d_o_mm')
    params.orf.d_o = params.orf.d_o_mm/1000;
end

req = {'Dp','d_w','D_m','n_turn','mu_ref','Lori','Lgap','Kd','Ebody','Gsh'};
if ~all(isfield(params,req)) || ~isfield(params,'orf') || ~isfield(params.orf,'d_o')
    return;
end

if ~isfield(params,'n_orf') && isfield(params,'orf') && isfield(params.orf,'n_orf')
    params.n_orf = params.orf.n_orf;
end
if ~isfield(params,'n_orf'), params.n_orf = 1; end
nd = 1;
if isfield(params,'nd') && isfinite(params.nd), nd = params.nd; end
if isfield(params,'n_parallel') && isfinite(params.n_parallel), nd = params.n_parallel; end
if isfield(params,'n_dampers_per_story')
    nds = params.n_dampers_per_story;
    if isnumeric(nds)
        if isscalar(nds), nd = nds; else, nd = max(1, round(max(nds(:)))); end
    end
end
nd = max(1, round(nd));

Ap = pi * params.Dp^2 / 4;
Ao_single = pi * params.orf.d_o^2 / 4;
Ao = params.n_orf * Ao_single;
Ap_eff = nd * Ap;
Ao_eff = nd * Ao;

rho_loc = util_getfield_default(params,'rho',850);
Lh = rho_loc * params.Lori / max(Ao^2, 1e-18);

k_h = params.Kd * Ap^2 / params.Lgap;
k_s = params.Ebody * Ap / params.Lgap;
k_hyd = 1 / (1/k_h + 1/k_s);
k_p = params.Gsh * params.d_w^4 / (8 * params.n_turn * params.D_m^3);
k_sd_simple = k_hyd + k_p;
k_sd_adv    = nd * (k_hyd + k_p);

c_lam0 = 12 * params.mu_ref * params.Lori * Ap^2 / (params.orf.d_o^4);

params.Ap = Ap;
params.Ao = Ao;
params.Ap_eff = Ap_eff;
params.Ao_eff = Ao_eff;
params.Lh = Lh;
params.k_p = k_p;
params.k_sd_simple = k_sd_simple;
params.k_sd_adv = k_sd_adv;
params.k_sd = k_sd_adv;
params.c_lam0 = c_lam0;
end

% UTIL_MAKE_ARIAS_WINDOW fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   t [-]: SI uyumlu giriş parametresi.
%   ag [-]: SI uyumlu giriş parametresi.
%   varargin [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   win [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function win = util_make_arias_window(t, ag, varargin)
%UTIL_MAKE_ARIAS_WINDOW Arias-intensity-based window generator.
p = inputParser;
p.addParameter('p1',0.05,@(x)isscalar(x) && x>=0 && x<=1);
p.addParameter('p2',0.95,@(x)isscalar(x) && x>=0 && x<=1);
p.addParameter('pad',0.5,@(x)isscalar(x) && x>=0);
p.parse(varargin{:});
p1 = p.Results.p1; p2 = p.Results.p2; pad = p.Results.pad;
IA = cumtrapz(t, ag.^2);
IA_tot = IA(end);
if ~(isfinite(IA_tot)) || IA_tot <= eps
    t_start = t(1); t_end = t(end);
    idx = true(size(t));
    win = struct('t5',t_start,'t95',t_end,'pad',0, ...
                 't_start',t_start,'t_end',t_end,'idx',idx, ...
                 'coverage',1.0,'flag_low_arias',true);
    return;
end
IA_norm = IA / IA_tot;
[t_unique, iu] = unique(IA_norm);
interp_t = t(iu);
t5  = interp1(t_unique, interp_t, p1, 'linear');
t95 = interp1(t_unique, interp_t, p2, 'linear');
dur = t95 - t5;
if dur < 5, pad = 0.25; end
t_start = max(t(1),  t5  - pad);
t_end   = min(t(end), t95 + pad);
idx = (t >= t_start) & (t <= t_end);
coverage = trapz(t(idx), ag(idx).^2) / IA_tot;
flag_low_arias = coverage < 0.90;
win = struct('t5',t5,'t95',t95,'pad',pad, ...
             't_start',t_start,'t_end',t_end,'idx',idx, ...
             'coverage',coverage,'flag_low_arias',flag_low_arias);
end

% UTIL_QUANTIZE_STEP fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   x [-]: SI uyumlu giriş parametresi.
%   step [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   y [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function y = util_quantize_step(x, step)
%UTIL_QUANTIZE_STEP Snap values to a uniform grid.
if nargin < 2 || isempty(step)
    y = x;
    return;
end
y = step * round(x ./ step);
end

% UTIL_GETFIELD_DEFAULT fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   S [-]: SI uyumlu giriş parametresi.
%   fname [-]: SI uyumlu giriş parametresi.
%   defaultVal [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   v [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function v = util_getfield_default(S, fname, defaultVal)
%UTIL_GETFIELD_DEFAULT Safe field accessor with defaults.
if ~isstruct(S) || ~isfield(S, fname) || isempty(S.(fname))
    v = defaultVal;
    return;
end
val = S.(fname);
if isnumeric(val)
    if isempty(val) || any(~isfinite(val(:)))
        v = defaultVal;
    else
        v = val;
    end
else
    v = val;
end
end

% UTIL_DEFAULT_PENALTY_OPTS fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   optsEval [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   lambda [-]: Hesaplanmış çıktı.
%   pwr [-]: Hesaplanmış çıktı.
%   W [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function [lambda,pwr,W] = util_default_penalty_opts(optsEval)
% util_default_penalty_opts — Ceza çekirdeğini tekleştirir.
% Varsayılanlar: λ=3 [-]: ölçek, pwr=2 [-]: üstel, W={dP=1,Qcap=1,cav=1,T=0.5,μ=0.5} [-]: ağırlıklar.
% Kullanıcı optsEval.penalty.(lambda|power|W) veya eski alanlar (penalty_scale,
%   penalty_power, penalty_weights) ile override edebilir.
    lambda = 3;
    pwr = 2;
    W = struct('dP',1,'Qcap',1,'cav',0.01,'T',0.5,'mu',0.5);
    if nargin < 1 || isempty(optsEval) || ~isstruct(optsEval)
        return;
    end

    if isfield(optsEval,'penalty') && isstruct(optsEval.penalty)
        pen_loc = optsEval.penalty;
        lambda = util_getfield_default(pen_loc,'lambda',lambda);
        pwr    = util_getfield_default(pen_loc,'power',pwr);
        if isfield(pen_loc,'W') && isstruct(pen_loc.W)
            W = merge_penalty_weights(W, pen_loc.W);
        end
    end

    if isfield(optsEval,'penalty_scale') && ~isempty(optsEval.penalty_scale)
        lambda = optsEval.penalty_scale;
    end
    if isfield(optsEval,'penalty_power') && ~isempty(optsEval.penalty_power)
        pwr = optsEval.penalty_power;
    end
    if isfield(optsEval,'penalty_weights') && isstruct(optsEval.penalty_weights)
        W = merge_penalty_weights(W, optsEval.penalty_weights);
    end
end

% MERGE_PENALTY_WEIGHTS fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   Wbase [-]: SI uyumlu giriş parametresi.
%   Woverride [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   Wout [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function Wout = merge_penalty_weights(Wbase, Woverride)
    Wout = Wbase;
    if ~isstruct(Woverride)
        return;
    end
    fn = fieldnames(Woverride);
    for k = 1:numel(fn)
        val_loc = Woverride.(fn{k});
        if isnumeric(val_loc) && isscalar(val_loc)
            Wout.(fn{k}) = val_loc;
        end
    end
end

% COMPUTE_THR fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   params [-]: SI uyumlu giriş parametresi.
%   thr_in [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   thr [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function thr = compute_thr(params, thr_in)
% compute_thr — QC eşiklerini tek kaynaktan üretir (override destekli).
% Hidrolik basınç üst sınırı Δp95_max [Pa], kapasite oranı Qcap95_max [-]:,
% kavitasyon yüzdesi cav_pct_max [-]:, son sıcaklık T_end_max [°C] ve son yağ
% viskozitesi μ_end_min [Pa·s] fiziksel kabullere dayanır. Kullanıcı thr_in
% alanları ile override yapabilir; eksik alanlar varsayılanla doldurulur.
    if nargin < 2 || ~isstruct(thr_in) || isempty(thr_in)
        thr_in = struct();
    end

    thermal_loc = util_getfield_default(params,'thermal', struct());
    cfg_loc     = util_getfield_default(params,'cfg', struct());
    num_loc     = util_getfield_default(cfg_loc,'num', struct());

    T0_C_loc    = util_getfield_default(params,'T0_C',25);
    dT_max_loc  = util_getfield_default(thermal_loc,'dT_max',80);
    mu_phys_loc = util_getfield_default(num_loc,'mu_min_phys',0.05);

    thr = thr_in;
    thr.dP95_max    = util_getfield_default(thr,'dP95_max',   40e6);  % [Pa] 40 MPa sınırı
    thr.Qcap95_max  = util_getfield_default(thr,'Qcap95_max', 0.90);  % [-]: kapasite oranı
    thr.cav_pct_max = util_getfield_default(thr,'cav_pct_max',0.00);  % [-]: kavitasyon yok
    thr.T_end_max   = util_getfield_default(thr,'T_end_max',  T0_C_loc + dT_max_loc); % [°C]
    thr.mu_end_min  = util_getfield_default(thr,'mu_end_min', max(mu_phys_loc, 0.6)); % [Pa·s]
end

% UTIL_INITIAL_POP_GRID fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   lb [-]: SI uyumlu giriş parametresi.
%   ub [-]: SI uyumlu giriş parametresi.
%   N [-]: SI uyumlu giriş parametresi.
%   steps [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   P [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function P = util_initial_pop_grid(lb, ub, N, steps)
%UTIL_INITIAL_POP_GRID Generate a lattice-aligned initial GA population.
d = numel(lb);
P = zeros(N, d);
for i = 1:d
    if ~isnan(steps(i))
        gs = steps(i);
        if gs <= 0
            gs = max(ub(i)-lb(i), eps)/10;
        end
        grid = lb(i):gs:ub(i);
        if numel(grid) < 2
            k = max(2, ceil((ub(i)-lb(i))/max(gs, eps)));
            grid = linspace(lb(i), ub(i), k);
        end
        idx = randi(numel(grid), [N,1]);
        vals = grid(idx);
    else
        vals = lb(i) + rand(N,1).*(ub(i)-lb(i));
    end
    P(:,i) = vals;
end
end

% UTIL_COMPUTE_CTH_EFFECTIVE fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   params [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   Cth [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function Cth = util_compute_Cth_effective(params)
%UTIL_COMPUTE_CTH_EFFECTIVE Effective combined thermal capacity.
nStories = size(params.M,1) - 1;
mask = params.story_mask(:);
if numel(mask)==1, mask = mask*ones(nStories,1); end
ndps = params.n_dampers_per_story(:);
if numel(ndps)==1, ndps = ndps*ones(nStories,1); end
multi = (mask .* ndps);
V_oil_per = params.resFactor * (params.Ap * (2*params.Lgap));
m_oil_tot = sum(multi) * (params.rho * V_oil_per);
m_steel_tot = params.steel_to_oil_mass_ratio * m_oil_tot;
Cth = max(m_oil_tot*params.cp_oil + m_steel_tot*params.cp_steel, eps);
end

% RESOLVE_THR_BASE fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   optsEval [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   thr [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function thr = resolve_thr_base(optsEval)
% resolve_thr_base — QC eşiklerini [Pa, -, °C, Pa·s] tek kaynaktan üretir.
% Öncelik sırası: optsEval.thr override, util_default_qc_thresholds() çıktısı,
% nihai fallback sabitleri (Qcap95_max [-]: 0.99, cav_pct_max [-]: 0.02,
% T_end_max [°C]: 75, μ_end_min [Pa·s]: 0.90, Δp95_max [Pa]: 30e6).
    if nargin < 1 || ~isstruct(optsEval)
        optsEval = struct();
    end

    thr = util_getfield_default(optsEval,'thr', struct());
    if ~isstruct(thr) || isempty(thr)
        thr = struct();
    end

    if isempty(fieldnames(thr))
        if exist('util_default_qc_thresholds','file') == 2
            thr_default = util_default_qc_thresholds();
            if isstruct(thr_default)
                thr = thr_default;
            end
        end
    end

    if ~isstruct(thr)
        thr = struct();
    end

    thr.Qcap95_max = util_getfield_default(thr,'Qcap95_max', 0.9);
    thr.cav_pct_max= util_getfield_default(thr,'cav_pct_max',0.00);
    thr.T_end_max  = util_getfield_default(thr,'T_end_max',  75);
    thr.mu_end_min = util_getfield_default(thr,'mu_end_min', 0.8);
    thr.dP95_max   = util_getfield_default(thr,'dP95_max',   40e6);
end
%% GA Sonrası Grafikler
% CIZ_GA_GRAFIKLERI fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   X [-]: SI uyumlu giriş parametresi.
%   scaled [-]: SI uyumlu giriş parametresi.
%   params [-]: SI uyumlu giriş parametresi.
%   optsEval [-]: SI uyumlu giriş parametresi.
% Çıktılar: Yok.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function ciz_ga_grafikleri(X, scaled, params, optsEval)
    if nargin < 1 || isempty(X) || isempty(scaled) || isempty(params)
        return;
    end
    if nargin < 4 || isempty(optsEval)
        optsEval = struct();
    end
    try
        best_design = X(1,:);
    catch
        return;
    end
    try
        best_params = decode_params_from_x(params, best_design);
        rec = scaled(1);
        assignin('base','plotx_design',best_design);
        assignin('base','plotx_params',best_params);
        assignin('base','plotx_rec',rec);
        sdof_like_dump_to_base(rec, best_params);
        plot_ga_sdof_comparisons(rec, best_params);
    catch ME
        warning('viskozmakale3:plot', 'GA grafik cizimi basarisiz: %s', ME.message);
    end
end

% PLOT_GA_SDOF_COMPARISONS fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   rec [-]: SI uyumlu giriş parametresi.
%   params [-]: SI uyumlu giriş parametresi.
% Çıktılar: Yok.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function plot_ga_sdof_comparisons(rec, params)
    if isempty(rec) || isempty(params)
        return;
    end

    t = rec.t(:);
    ag = rec.ag(:);
    n = size(params.M, 1);
    iTop = n;

    win = util_make_arias_window(t, ag);

    [x_linear, a_rel_linear] = ga_util_solve_linear_mck(t, ag, params.M, params.C0, params.K);

    mu_ref_eff = params.mu_ref;
    c_lam0_eff = params.c_lam0;
    if isfield(params, 'cfg') && isstruct(params.cfg) && isfield(params.cfg, 'on') ...
            && isstruct(params.cfg.on) && isfield(params.cfg.on, 'mu_floor') && params.cfg.on.mu_floor
        mu_min_phys = NaN;
        if isfield(params.cfg, 'num') && isstruct(params.cfg.num) && ...
                isfield(params.cfg.num, 'mu_min_phys') && isfinite(params.cfg.num.mu_min_phys)
            mu_min_phys = params.cfg.num.mu_min_phys;
        end
        if isfinite(mu_min_phys) && mu_ref_eff < mu_min_phys
            mu_ref_eff = mu_min_phys;
            scale_mu = mu_ref_eff / max(params.mu_ref, eps);
            c_lam0_eff = params.c_lam0 * scale_mu;
        end
    end

    Tinit = params.T_ref_C;

    [x_damper, a_rel_damper] = mck_with_damper(t, ag, params.M, params.C0, params.K, ...
        params.k_sd, c_lam0_eff, params.Lori, params.orf, params.rho, ...
        params.Ap, params.Ao, params.Qcap_big, mu_ref_eff, params.thermal, ...
        Tinit, params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, ...
        params.Lgap, params.cp_oil, params.cp_steel, ...
        params.steel_to_oil_mass_ratio, params.story_mask, ...
        params.n_dampers_per_story, params.resFactor, params.cfg);

    x_top_linear = x_linear(:, iTop);
    x_top_damper = x_damper(:, iTop);
    a_abs_linear = a_rel_linear(:, iTop) + ag;
    a_abs_damper = a_rel_damper(:, iTop) + ag;

    T1 = ga_util_compute_T1(params.M, params.K);

    figure('Name', sprintf('Ust Kat Yer Degistirme - %s', rec.name), 'Color', 'w');
    plot(t, x_top_linear, 'k', 'LineWidth', 1.4); hold on;
    plot(t, x_top_damper, 'r', 'LineWidth', 1.0);
    yl = ylim;
    plot([win.t5 win.t5], yl, 'k--', 'HandleVisibility', 'off');
    plot([win.t95 win.t95], yl, 'k--', 'HandleVisibility', 'off');
    grid on;
    xlabel('Zaman [s]');
    ylabel(sprintf('x_{%d}(t) [m]', iTop));
    title(sprintf('Ust Kat Yer Degistirme | T1 = %.3f s | Arias [%.3f, %.3f] s | %s', ...
        T1, win.t5, win.t95, rec.name), 'Interpreter', 'none');
    legend('Dampersiz', 'Damperli', 'Location', 'best');

    figure('Name', sprintf('Ust Kat Mutlak Ivme - %s', rec.name), 'Color', 'w');
    plot(t, a_abs_linear, 'k', 'LineWidth', 1.4); hold on;
    plot(t, a_abs_damper, 'r', 'LineWidth', 1.0);
    yl = ylim;
    plot([win.t5 win.t5], yl, 'k--', 'HandleVisibility', 'off');
    plot([win.t95 win.t95], yl, 'k--', 'HandleVisibility', 'off');
    grid on;
    xlabel('Zaman [s]');
    ylabel('a_{abs}(t) [m/s^2]');
    title(sprintf('Ust Kat Mutlak Ivme | T1 = %.3f s | Arias [%.3f, %.3f] s | %s', ...
        T1, win.t5, win.t95, rec.name), 'Interpreter', 'none');
    legend('Dampersiz', 'Damperli', 'Location', 'best');

    drift_linear = x_linear(:, 2:end) - x_linear(:, 1:end-1);
    drift_damper = x_damper(:, 2:end) - x_damper(:, 1:end-1);
    idx = win.idx;
    IDR_linear = max(abs(drift_linear(idx, :)), [], 1) ./ params.story_height;
    IDR_damper = max(abs(drift_damper(idx, :)), [], 1) ./ params.story_height;
    story_ids = 1:(n-1);

    figure('Name', sprintf('Maksimum IDR - %s', rec.name), 'Color', 'w');
    plot(story_ids, IDR_linear, 'k-o', 'LineWidth', 1.4); hold on;
    plot(story_ids, IDR_damper, 'r-d', 'LineWidth', 1.0);
    grid on;
    xlabel('Kat');
    ylabel('Maksimum IDR [\Delta x / h]');
    title(sprintf('Katlar Arasi Maksimum Otelenme Orani | %s', rec.name), 'Interpreter', 'none');
    legend('Dampersiz', 'Damperli', 'Location', 'best');
end

% GA_UTIL_SOLVE_LINEAR_MCK fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   t [-]: SI uyumlu giriş parametresi.
%   ag [-]: SI uyumlu giriş parametresi.
%   M [-]: SI uyumlu giriş parametresi.
%   C [-]: SI uyumlu giriş parametresi.
%   K [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   x [-]: Hesaplanmış çıktı.
%   a_rel [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function [x, a_rel] = ga_util_solve_linear_mck(t, ag, M, C, K)
    nloc = size(M, 1);
    r = ones(nloc, 1);
    agf = griddedInterpolant(t, ag, 'linear', 'nearest');
    odefun = @(tt, z) [z(nloc+1:end); ...
        M \ (-C * z(nloc+1:end) - K * z(1:nloc) - M * r * agf(tt))];
    z0 = zeros(2 * nloc, 1);
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
    sol = ode15s(odefun, [t(1) t(end)], z0, opts);
    z = deval(sol, t).';
    x = z(:, 1:nloc);
    v = z(:, nloc+1:end);
    a_rel = (-(M \ (C * v.' + K * x.')).' - ag(:) .* r.');
end

% GA_UTIL_COMPUTE_T1 fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   M [-]: SI uyumlu giriş parametresi.
%   K [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   T1 [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function T1 = ga_util_compute_T1(M, K)
    [~, D] = eig(K, M);
    w = sqrt(sort(diag(D), 'ascend'));
    T1 = 2 * pi / w(1);
end

% SDOF_LIKE_DUMP_TO_BASE fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   rec [-]: SI uyumlu giriş parametresi.
%   params [-]: SI uyumlu giriş parametresi.
% Çıktılar: Yok.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function sdof_like_dump_to_base(rec, params)
    t  = rec.t(:);
    ag = rec.ag(:);
    n  = size(params.M,1);
    iTop = n;

    win = util_make_arias_window(t, ag);
    t5 = win.t5; t95 = win.t95;

    [x0, a_rel0] = util_solve_linear_mck(t, ag, params.M, params.C0, params.K);

    mu_ref_eff = params.mu_ref;
    c_lam0_eff = params.c_lam0;
    if isfield(params,'cfg') && isstruct(params.cfg) && isfield(params.cfg,'on') ...
            && isfield(params.cfg.on,'mu_floor') && params.cfg.on.mu_floor
        mu_min_phys = NaN;
        if isfield(params,'cfg') && isfield(params.cfg,'num') && isfield(params.cfg.num,'mu_min_phys')
            mu_min_phys = params.cfg.num.mu_min_phys;
        end
        if isfinite(mu_min_phys) && (mu_ref_eff < mu_min_phys)
            mu_ref_eff = mu_min_phys;
            scl = mu_ref_eff / max(params.mu_ref, eps);
            c_lam0_eff = params.c_lam0 * scl;
        end
    end
    Tinit = params.T_ref_C;

    [x_d, a_d, diag_d] = mck_with_damper(t, ag, params.M, params.C0, params.K, ...
        params.k_sd, c_lam0_eff, params.Lori, params.orf, params.rho, ...
        params.Ap, params.Ao, params.Qcap_big, mu_ref_eff, params.thermal, ...
        Tinit, params.T_ref_C, params.b_mu, params.c_lam_min, params.c_lam_cap, ...
        params.Lgap, params.cp_oil, params.cp_steel, ...
        params.steel_to_oil_mass_ratio, params.story_mask, ...
        params.n_dampers_per_story, params.resFactor, params.cfg);

    x10_0 = x0(:,iTop);    x10_d = x_d(:,iTop);
    a10_0 = a_rel0(:,iTop) + ag;
    a10_d = a_d(:,iTop)    + ag;

    drift0 = x0(:,2:end) - x0(:,1:end-1);
    drift_d = x_d(:,2:end) - x_d(:,1:end-1);
    idx = win.idx;
    story_height = params.story_height;
    IDR0 = max(abs(drift0(idx,:)))./story_height;
    IDR_d = max(abs(drift_d(idx,:)))./story_height;

    Ns = n - 1; mask = params.story_mask(:); if numel(mask)==1, mask = mask*ones(Ns,1); end
    ndps = params.n_dampers_per_story(:);    if numel(ndps)==1, ndps = ndps*ones(Ns,1); end
    multi = (mask .* ndps);

    Kadd = zeros(n); C_add = zeros(n);
    for i = 1:Ns
        ii = [i,i+1];
        k_eq = params.k_sd * multi(i);
        c_eq = diag_d.c_lam * multi(i);
        Kadd(ii,ii) = Kadd(ii,ii) + k_eq * [ 1 -1; -1 1 ];
        C_add(ii,ii)= C_add(ii,ii) + c_eq * [ 1 -1; -1 1 ];
    end
    K_tot = params.K + Kadd;
    C_d   = params.C0 + C_add;

    [V,D] = eig(K_tot, params.M); [w2,ord] = sort(diag(D),'ascend');
    phi1 = V(:,ord(1)); w1 = sqrt(w2(1));
    normM = phi1.' * params.M * phi1;
    T1 = util_compute_T1(params.M, params.K);
    zeta0 = (phi1.' * params.C0 * phi1) / (2*w1*normM);
    zeta_d = (phi1.' * C_d        * phi1) / (2*w1*normM);

    x10_max_0            = max(abs(x10_0(idx)));
    x10_max_damperli     = max(abs(x10_d(idx)));
    a10abs_max_0         = max(abs(a10_0(idx)));
    a10abs_max_damperli  = max(abs(a10_d(idx)));
    IDR_max_0            = max(IDR0);
    IDR_max_damperli     = max(IDR_d);
    story_ids            = 1:Ns;

    assignin('base','t',t);           assignin('base','ag',ag);
    assignin('base','win',win);       assignin('base','t5',t5);      assignin('base','t95',t95);
    assignin('base','x0',x0);         assignin('base','a_rel0',a_rel0);
    assignin('base','x_d',x_d);       assignin('base','a_d',a_d);    assignin('base','diag_d',diag_d);
    assignin('base','x10_0',x10_0);   assignin('base','x10_d',x10_d);
    assignin('base','a10_0',a10_0);   assignin('base','a10_d',a10_d);
    assignin('base','drift0',drift0); assignin('base','drift_d',drift_d);
    assignin('base','IDR0',IDR0);     assignin('base','IDR_d',IDR_d);
    assignin('base','story_ids',story_ids);
    assignin('base','K_tot',K_tot);   assignin('base','C_d',C_d);
    assignin('base','phi1',phi1);     assignin('base','w1',w1);      assignin('base','T1',T1);
    assignin('base','zeta0',zeta0);   assignin('base','zeta_d',zeta_d);
    assignin('base','x10_max_0',x10_max_0);
    assignin('base','x10_max_damperli',x10_max_damperli);
    assignin('base','a10abs_max_0',a10abs_max_0);
    assignin('base','a10abs_max_damperli',a10abs_max_damperli);
    assignin('base','IDR_max_0',IDR_max_0);
    assignin('base','IDR_max_damperli',IDR_max_damperli);
end

% UTIL_SOLVE_LINEAR_MCK fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   varargin [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   x [-]: Hesaplanmış çıktı.
%   a_rel [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function [x, a_rel] = util_solve_linear_mck(varargin)
    [x, a_rel] = ga_util_solve_linear_mck(varargin{:});
end

% UTIL_COMPUTE_T1 fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   varargin [-]: SI uyumlu giriş parametresi.
% Çıktılar:
%   T1 [-]: Hesaplanmış çıktı.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function T1 = util_compute_T1(varargin)
    T1 = ga_util_compute_T1(varargin{:});
end


% RAPORLA_GA_SONUCLARI fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler:
%   X [-]: SI uyumlu giriş parametresi.
%   F [-]: SI uyumlu giriş parametresi.
%   scaled [-]: SI uyumlu giriş parametresi.
%   params [-]: SI uyumlu giriş parametresi.
%   optsEval [-]: SI uyumlu giriş parametresi.
%   gaout [-]: SI uyumlu giriş parametresi.
% Çıktılar: Yok.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function raporla_ga_sonuclari(X, F, scaled, params, optsEval, gaout)
    if nargin < 1 || isempty(X)
        if nargin < 3 || isempty(scaled) || isempty(params)
            [scaled, params] = prepare_inputs();
        end
        if nargin < 5 || isempty(optsEval)
            optsEval = struct();
        end
        fprintf('Kayitli Pareto sonucu bulunamadi. GA optimizasyonu baslatiliyor...\n');
        [X, F, gaout] = run_ga_driver(scaled, params, optsEval, struct());
    else
        if nargin < 3 || isempty(scaled) || isempty(params)
            [scaled, params] = prepare_inputs();
        end
        if nargin < 2 || isempty(F)
            F = [];
        end
        if nargin < 5 || isempty(optsEval)
            optsEval = struct();
        end
        if nargin < 6 || isempty(gaout)
            gaout = struct();
        end
    end

    if isempty(X)
        warning('viskozmakale3:rapor', 'GA cikisinda gecerli tasarim yok.');
        return;
    end

    fprintf('Pareto on cephesindeki tasarim sayisi: %d\n', size(X, 1));

    best_design = X(1, :);
    [fitness_best, ~, details_best] = eval_design_fast(best_design, scaled, params, optsEval);
    fprintf('En iyi aday hedefleri (pen, f1, f2): %.4g / %.4g / %.4g\n', fitness_best(1), fitness_best(2), fitness_best(3));

    fprintf('Karar vektoru degerleri:\n');
    fprintf('  d_o_mm = %.2f, n_orf = %d, Cd0 = %.2f, CdInf = %.2f, p_exp = %.2f\n', best_design(1), round(best_design(2)), best_design(3), best_design(4), best_design(5));
    fprintf('  Lori_mm = %.1f, hA_W_perK = %.1f, Dp_mm = %.1f, d_w_mm = %.1f, D_m_mm = %.1f\n', best_design(6), best_design(7), best_design(8), best_design(9), best_design(10));
    fprintf('  n_turn = %d, mu_ref = %.2f\n', round(best_design(11)), best_design(12));

    best_params = decode_params_from_x(params, best_design);
    if isfield(best_params, 'thermal') && isfield(best_params.thermal, 'hA_W_perK')
        fprintf('Termal kapasite hA_W_perK (W/K): %.2f\n', best_params.thermal.hA_W_perK);
    end

    if ~isempty(F)
        fprintf('Kayitli uygunluk degerleri (ilk birey): %.4g / %.4g / %.4g\n', F(1, 1), F(1, 2), F(1, 3));
    end

    if exist('details_best', 'var') && isstruct(details_best) && isfield(details_best, 'table')
        tablo = details_best.table;
        alanlar = {'PFA','IDR','dP95','Qcap95','cav_pct','T_end','mu_end'};
        fprintf('Ozet metrikler (ortalama):\n');
        for i = 1:numel(alanlar)
            ad = alanlar{i};
            if isfield(tablo, ad)
                val = mean(enforce_finite(tablo.(ad)));
                fprintf('  %s: %.4g\n', ad, val);
            end
        end
    end

    if isstruct(gaout) && isfield(gaout, 'exitflag')
        fprintf('GA cikis kodu: %d\n', gaout.exitflag);
    end
end

% RAPORLA_KAYITLI_GA_SONUCLARI fonksiyonu için otomatik dokümantasyon bloğu.
% Amaç ve Kapsam: Bu yardımcı fonksiyon, viskoz damper GA akışındaki ilgili hesap adımını yürütür.
% Girdiler: Yok.
% Çıktılar: Yok.
% Varsayımlar ve Sınırlar: - Girdi ve çıktılar tutarlı boyutlara sahip olmalıdır.
% Ölçüler/Birimler: SI tabanlı varsayılır; özgül birimler kullanım bağlamında belirlenir.
% Yöntem Özeti: Fonksiyonun iç mantığı mevcut yardımcılarla tanımlıdır.
function raporla_kayitli_ga_sonuclari()
    ga_dirs = dir('out/ga_*');
    if isempty(ga_dirs)
        raporla_ga_sonuclari([], [], [], [], struct(), struct());
        return;
    end

    [~, idx] = max([ga_dirs.datenum]);
    latest_dir = ga_dirs(idx).name;
    fprintf('Kayitli GA sonucu yukleniyor: %s\n', latest_dir);

    ga_data = load(fullfile('out', latest_dir, 'ga_front.mat'));

    if ~isfield(ga_data, 'X') || isempty(ga_data.X)
        raporla_ga_sonuclari([], [], [], [], struct(), struct());
        return;
    end

    [scaled, params] = prepare_inputs();

    storedF = [];
    if isfield(ga_data, 'F')
        storedF = ga_data.F;
    end

    raporla_ga_sonuclari(ga_data.X, storedF, scaled, params, struct(), struct());
end
