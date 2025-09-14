% Bu dosya, farklı statik yardımcı fonksiyonları barındıran Utils sınıfını tanımlar.
classdef Utils
%UTILS Küçük yardımcı fonksiyonları statik yöntemler olarak toplayan sınıf.
    methods(Static)
        %% Yumuşak Minimum
        function y = softmin(a,b,epsm)
            % İki değerin yumuşak minimumunu hesaplar.
            % Örnek kullanım: y = Utils.softmin(3,5,0.2);
            y = 0.5*(a + b - sqrt((a - b).^2 + epsm.^2));
        end

        %% Softmin epsilon (�l�ekli)
        function epsm = softmin_eps(cfg)
            % SOFTMIN_EPS  dP yumu�atma epsilonu (�l�ekli) d�nd�r�r.
            %  epsm = c_eps * num.dP_cap; c_eps varsay�lan 0.03.
            try
                c_eps = Utils.getfield_default(cfg.num,'softmin_ceps',0.03);
                dPcap = Utils.getfield_default(cfg.num,'dP_cap',NaN);
                if isfinite(dPcap) && dPcap>0
                    epsm = max(1e3, c_eps * dPcap);
                else
                    epsm = 1e5; % emniyetli varsay�lan
                end
            catch
                epsm = 1e5;
            end
        end

        %% Orantılı Pencere Ağırlığı
        function w = pf_weight(t, cfg)
            % Basınç kuvveti için orantılı pencere ağırlığı hesaplar.
            % Örnek kullanım: w = Utils.pf_weight(t, cfg);
            % Tek-kaynak PF ramp: compat_simple (eski) vs softplus (ileri)
            if ~isstruct(cfg), cfg = struct(); end
            if ~isfield(cfg,'on') || ~isstruct(cfg.on), cfg.on = struct(); end
            if ~isfield(cfg.on,'pressure_force'), cfg.on.pressure_force = true; end
            if ~isfield(cfg,'PF') || ~isstruct(cfg.PF), cfg.PF = struct(); end
            if ~isfield(cfg.PF,'t_on'), cfg.PF.t_on = 0; end
            if ~isfield(cfg.PF,'tau'),  cfg.PF.tau  = 1.0; end
            if ~isfield(cfg,'compat_simple'), cfg.compat_simple = true; end

            t = double(t); tau_floor = 1e-6;
            if cfg.compat_simple
                dt  = max(t - cfg.PF.t_on, 0);
                tau = max(cfg.PF.tau, tau_floor);
                w_local = 1 - exp(-dt ./ tau);
            else
                k = Utils.getfield_default(cfg.PF,'k', 0.01); k = max(k, tau_floor);
                sp_dt  = (log1p(exp(-abs((t - cfg.PF.t_on)./k))) + max((t - cfg.PF.t_on)./k, 0));
                dt  = sp_dt .* k;
                sp_tau = (log1p(exp(-abs((cfg.PF.tau - tau_floor)./k))) + max((cfg.PF.tau - tau_floor)./k, 0));
                tau = sp_tau .* k + tau_floor;
                w_local = 1 - exp(-dt ./ tau);
            end
            w = cfg.on.pressure_force .* w_local;
        end

        %% Damper Sabitlerini Güncelle
        function params = recompute_damper_params(params)
            %RECOMPUTE_DAMPER_PARAMS Türetilmiş damper sabitlerini günceller.
            %   PARAMS = RECOMPUTE_DAMPER_PARAMS(PARAMS) yapısı içindeki
            %   temel geometrik ve malzeme parametrelerine (Dp, d_w, D_m,
            %   n_turn, mu_ref vb.) göre Ap, k_p, k_sd ve c_lam0 gibi
            %   türetilmiş sabitleri yeniden hesaplar. Eksik alanlar
            %   bulunduğunda mevcut değerler korunur.

            if ~isstruct(params), return; end

            % mm cinsinden verilen değerleri metreye çevir
            if isfield(params,'Dp_mm'),    params.Dp    = params.Dp_mm/1000; end
            if isfield(params,'d_w_mm'),   params.d_w   = params.d_w_mm/1000; end
            if isfield(params,'D_m_mm'),   params.D_m   = params.D_m_mm/1000; end
            if isfield(params,'Lori_mm'),  params.Lori  = params.Lori_mm/1000; end
            if isfield(params,'orf') && isfield(params.orf,'d_o_mm')
                params.orf.d_o = params.orf.d_o_mm/1000;
            end

            req = {'Dp','d_w','D_m','n_turn','mu_ref','Lori','Lgap','Kd','Ebody','Gsh'};
            if ~all(isfield(params,req)) || ~isfield(params,'orf') || ~isfield(params.orf,'d_o')
                return; % eksik alanlar varsa hesaplama yapma
            end

            % n_orf ve nd (paralel) güvenli varsayılanlar
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

            % Alanlar (tek damper ve efektif)
            Ap = pi * params.Dp^2 / 4;
            Ao_single = pi * params.orf.d_o^2 / 4;
            Ao = params.n_orf * Ao_single;
            Ap_eff = nd * Ap;
            Ao_eff = nd * Ao;

            % Hat atalet Lh (tek damper)
            rho_loc = Utils.getfield_default(params,'rho',850);
            Lh = rho_loc * params.Lori / max(Ao^2, 1e-18);

            % Rijitlikler (tek damper)
            k_h = params.Kd * Ap^2 / params.Lgap;
            k_s = params.Ebody * Ap / params.Lgap;
            k_hyd = 1 / (1/k_h + 1/k_s);
            k_p = params.Gsh * params.d_w^4 / (8 * params.n_turn * params.D_m^3);
            k_sd_simple = k_hyd + k_p;      % tek damper
            k_sd_adv    = nd * (k_hyd + k_p);% paralel nd

            % Laminer sabit (tek damper referansı)
            c_lam0 = 12 * params.mu_ref * params.Lori * Ap^2 / (params.orf.d_o^4);

            % Çıkışlar (geriye uyumlu alan adlarıyla)
            params.Ap = Ap;
            params.Ao = Ao; params.A_o = Ao;
            params.Ap_eff = Ap_eff; params.Ao_eff = Ao_eff;
            params.Lh = Lh;
            params.k_p = k_p;
            params.k_sd_simple = k_sd_simple;
            params.k_sd_adv = k_sd_adv;
            % Adım 2 öncesi: k_sd paralel etkili (nd içselleştirilmiş) seçilir
            params.k_sd = k_sd_adv;

            params.c_lam0 = c_lam0;
        end

        %% Lineer MCK Çözümü
        function [x,a] = lin_MCK(t,ag,M,C,K)
            % Lineer MCK sistemi için yer hareketi altındaki tepkiyi çözer.
            % Örnek kullanım: [x,a] = Utils.lin_MCK(t, ag, M, C, K);
            n = size(M,1); r = ones(n,1);
            agf = griddedInterpolant(t,ag,'linear','nearest');
            odef = @(tt,z)[ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - M*r*agf(tt) ) ];
            z0 = zeros(2*n,1);
            opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
            sol = ode15s(odef,[t(1) t(end)],z0,opts);
            z = deval(sol,t).';
            x = z(:,1:n);
            a = ( -(M\(C*z(:,n+1:end).' + K*z(:,1:n).')).' - ag.*r.' );
        end

        %% Arias Penceresi Oluşturma
        function win = make_arias_window(t, ag, varargin)
            % Arias yoğunluğu tabanlı pencere oluşturur.
            % Örnek kullanım: win = Utils.make_arias_window(t, ag);
            p = inputParser;
            p.addParameter('p1',0.05,@(x)isscalar(x) && x>=0 && x<=1);
            p.addParameter('p2',0.95,@(x)isscalar(x) && x>=0 && x<=1);
            p.addParameter('pad',0.5,@(x)isscalar(x) && x>=0);
            p.parse(varargin{:});
            p1 = p.Results.p1; p2 = p.Results.p2; pad = p.Results.pad;
            IA = cumtrapz(t, ag.^2);
            IA_tot = IA(end);
            % Bozulmuş veya çok düşük enerjili kayıtlara karşı koruma
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

        %% Adımsal Nicemleme
        function y = quantize_step(x, step)
            % Verilen adım büyüklüğüne göre x değerini nicemler.
            % Örnek kullanım: y = Utils.quantize_step(3.7, 0.5);
            if nargin < 2 || isempty(step)
                y = x; return;
            end
            y = step * round(x ./ step);
        end

        %% Varsayılan Alan Değeri
        function v = getfield_default(S, fname, defaultVal)
            % Yapı alanı mevcut değilse varsayılan değeri döndürür.
            % Örnek kullanım: v = Utils.getfield_default(S,'a',0);
            if ~isstruct(S) || ~isfield(S, fname) || isempty(S.(fname))
                v = defaultVal; return;
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

        %% Alan Mevcutsa Atama
        function arr = assign_if_field(S, fname, arr, idx)
            % Belirtilen alan mevcutsa değeri hedef dizinin idx konumuna atar.
            % Örnek kullanım: Q = Utils.assign_if_field(m,'Q_q95',Q,k);
            if isfield(S, fname)
                arr(idx) = S.(fname);
            end
        end

        %% JSON Yazma
        function writejson(data, filename)
            % Verilen veriyi JSON dosyasına yazar.
            % Örnek kullanım: Utils.writejson(data,'cikti.json');
            try
                txt = jsonencode(data);
                fid = fopen(filename,'w');
                if fid~=-1
                    fwrite(fid, txt);
                    fclose(fid);
                end
            catch
            end
        end

        %% İsim Temizleme
        function s2 = sanitize_name(s)
            % Dosya veya alan isimlerindeki geçersiz karakterleri temizler.
            % Örnek kullanım: s2 = Utils.sanitize_name('örnek?*ad');
            if ~ischar(s) && ~isstring(s)
                s = char(s);
            end
            s2 = regexprep(char(s),'[^a-zA-Z0-9_\- ]','_');
        end

        %% Zaman Serisi Aşağı Örnekleme
        function ts_ds = downsample_ts(ts, ds)
            % Zaman serisi alanlarını verilen faktörle seyrekleştirir.
            % Örnek kullanım: ts_ds = Utils.downsample_ts(ts, 5);
            if nargin < 2 || isempty(ds), ds = 5; end
            fns = fieldnames(ts);
            ts_ds = struct();
            for i = 1:numel(fns)
                val = ts.(fns{i});
                if isnumeric(val) && ~isscalar(val)
                    idx = 1:ds:size(val,1);
                    ts_ds.(fns{i}) = val(idx,:,:);
                else
                    ts_ds.(fns{i}) = val;
                end
            end
        end

        %% Üçlü Operatör
        function s = tern(c,a,b)
            % Mantıksal koşula göre iki değerden birini seçer.
            % Örnek kullanım: s = Utils.tern(x>0, 1, -1);
            if c, s=a; else, s=b; end
        end

        %% Güvenli Çalıştırma ve Uyarı
        function varargout = try_warn(funcHandle, msgPrefix)
            %TRY_WARN Hata durumunda uyarı vererek fonksiyonu çalıştırır.
            %   varargout = Utils.try_warn(@() func(), 'prefix')
            try
                [varargout{1:nargout}] = funcHandle();
            catch ME
                warning('%s: %s', msgPrefix, ME.message);
                for k = 1:nargout
                    varargout{k} = [];
                end
            end
        end

        %% Başlangıç Popülasyon Izgarası
        function P = initial_pop_grid(lb, ub, N, steps)
            %INITIAL_POP_GRID Değişken ızgaralarına hizalanmış bir başlangıç popülasyonu oluşturur.
            %   P = Utils.initial_pop_grid(lb, ub, N, steps) fonksiyonu NxD boyutlu bir matris
            %   döndürür. Her i sütunu adım aralığı steps(i) olan bir ızgaradan örneklenir
            %   (steps(i) NaN ise, [lb(i), ub(i)] aralığında uniform örnekleme yapılır).
            % Örnek kullanım: P = Utils.initial_pop_grid([0 0],[1 1],5,[0.1 NaN]);
            d = numel(lb);
            P = zeros(N, d);
            for i = 1:d
                if ~isnan(steps(i))
                    if steps(i) <= 0
                        gs = max(ub(i)-lb(i), eps)/10;
                    else
                        gs = steps(i);
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

        %% Varsayılan QC Eşikleri
        function thr = default_qc_thresholds(optsThr)
            % Kalite kontrolü için varsayılan eşik değerlerini döndürür ve
            % eksik veya boş alanları varsayılanlarla doldurur.
            % Örnek kullanım: thr = Utils.default_qc_thresholds(struct('dP95_max',40e6));

            if nargin < 1 || isempty(optsThr)
                optsThr = struct();
            end

            thr = struct();
            thr.dP95_max   = Utils.getfield_default(optsThr,'dP95_max',50e6);
            thr.Qcap95_max = Utils.getfield_default(optsThr,'Qcap95_max',0.5);
            thr.cav_pct_max= Utils.getfield_default(optsThr,'cav_pct_max',0);
            thr.T_end_max  = Utils.getfield_default(optsThr,'T_end_max',75);
            thr.mu_end_min = Utils.getfield_default(optsThr,'mu_end_min',0.5);

            % Ek alanları koru
            extra = setdiff(fieldnames(optsThr), fieldnames(thr));
            for ii = 1:numel(extra)
                thr.(extra{ii}) = optsThr.(extra{ii});
            end
        end
    end
end
