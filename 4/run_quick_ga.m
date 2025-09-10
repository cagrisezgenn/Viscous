function [X,F,gaout,outdir] = run_quick_ga()
% RUN_QUICK_GA Hızlı başlangıç: veriyi ölçekle, GA'yı önerilen ayarlarla çalıştır.

    % Yol hazırlığı
    try, setup; catch, end

    % Temel parametreleri yükle ve T1 hesapla
    parametreler; %#ok<NODEF>

    % Kayıtları yükle ve ölçekle (band IM)
    L = struct('IM_mode','band','band_fac',[0.8 1.2],'s_bounds',[0.5 2.0]);
    [~, scaled, meta] = load_ground_motions(T1, L); %#ok<NASGU>

    % Parametre yapısını oluştur
    params = struct();
    params.M = M; params.C0 = C0; params.K = K;
    params.k_sd = k_sd; params.c_lam0 = c_lam0; params.Lori = Lori;
    params.orf = orf; params.rho = rho; params.Ap = Ap; params.A_o = A_o;
    params.Qcap_big = Qcap_big; params.mu_ref = mu_ref;
    params.thermal = thermal; params.T0_C = T0_C; params.T_ref_C = T_ref_C; params.b_mu = b_mu;
    params.c_lam_min = c_lam_min; params.c_lam_cap = c_lam_cap; params.Lgap = Lgap;
    params.cp_oil = cp_oil; params.cp_steel = cp_steel; params.steel_to_oil_mass_ratio = steel_to_oil_mass_ratio;
    params.n_dampers_per_story = n_dampers_per_story; params.toggle_gain = toggle_gain; params.story_mask = story_mask;
    params.story_height = story_height;  % gerekli saha: compute_metrics_windowed için
    params.resFactor = resFactor; params.cfg = cfg; params.n_orf = n_orf;
    params.Dp = Dp; params.d_w = d_w; params.D_m = D_m; params.n_turn = n_turn;
    params.Kd = Kd; params.Ebody = Ebody; params.Gsh = Gsh;

    % Değerlendirme ve GA opsiyonları (önerilen)
    optsEval = struct();
    optsEval.mu_factors = [0.70 1.00 1.30];
    optsEval.mu_weights = [0.3 0.5 0.2];
    optsEval.penalty_scale = 5;
    optsEval.penalty_power = 1.0;
    optsEval.penalty_weights = struct('dP',1,'Qcap',1,'cav',3,'T',1,'mu',0.5);
    optsEval.hardkill = struct('dP',1.5,'Qcap',1.2,'cav_abs',0.02);

    optsGA = struct();
    optsGA.PopulationSize  = 120;
    optsGA.MaxGenerations  = 80;
    optsGA.CrossoverFraction = 0.8;
    optsGA.ParetoFraction  = 0.50;
    optsGA.UseParallel     = true;

    % Çalıştır
    [X,F,gaout] = run_ga_driver(scaled, params, optsEval, optsGA);

    % Son çıkan klasörü yakala
    dd = dir(fullfile('out','ga_*'));
    if ~isempty(dd)
        [~,ix] = max([dd.datenum]);
        outdir = fullfile(dd(ix).folder, dd(ix).name);
    else
        outdir = '';
    end

    fprintf('GA tamamlandı. Sonuç klasörü: %s\n', outdir);
end
