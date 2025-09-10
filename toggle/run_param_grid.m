%% Parametre ızgarası üzerinde analiz çalıştır
% Bu betik parametreler.m içindeki tarama dizilerini alır ve
% tüm kombinasyonlar için run_one_record_windowed fonksiyonunu çağırır.

% 1) Parametreleri yükle
parametreler; %#ok<*NASGU>

% Tarama yapılacak parametre dizileri
ndps_vals = n_dampers_per_story(:);
gain_vals = toggle_gain(:);

% 2) Tüm kombinasyonları üret
[NDPS, GAIN] = ndgrid(ndps_vals, gain_vals);
comb = [NDPS(:), GAIN(:)];

% 3) Yer hareketi kaydını yükle (hızlı olması için ilk kayıt)
[~, scaled] = load_ground_motions(T1);
rec = scaled(1);

% 4) Temel parametre yapısı (tarama dışındaki parametreler)
base_params = struct('M',M,'C0',C0,'K',K,'k_sd',k_sd,'c_lam0',c_lam0, ...
    'orf',orf,'rho',rho,'Ap',Ap,'A_o',A_o,'Qcap_big',Qcap_big,'mu_ref',mu_ref, ...
    'thermal',thermal,'T0_C',T0_C,'T_ref_C',T_ref_C,'b_mu',b_mu, ...
    'c_lam_min',c_lam_min,'c_lam_cap',c_lam_cap,'Lgap',Lgap, ...
    'cp_oil',cp_oil,'cp_steel',cp_steel,'steel_to_oil_mass_ratio',steel_to_oil_mass_ratio, ...
    'story_mask',story_mask,'resFactor',resFactor,'cfg',cfg,'story_height',story_height);

% Çıktı dosyası
outdir = 'out'; if ~exist(outdir,'dir'), mkdir(outdir); end
outfile = fullfile(outdir, 'param_grid_results.csv');
if exist(outfile,'file'), delete(outfile); end
all_results = table();

% 5) Parametre kombinasyonları üzerinde döngü
opts = struct('use_orifice', true, 'use_thermal', true);
for i = 1:size(comb,1)
    params = base_params;
    params.n_dampers_per_story = comb(i,1);
    params.toggle_gain = comb(i,2);

    out = run_one_record_windowed(rec, [], params, opts);
    row = table(params.n_dampers_per_story, params.toggle_gain, ...
        out.metr.PFA_top, out.metr.IDR_max, out.metr.dP_orf_q95, ...
        'VariableNames', {'n_dampers_per_story','toggle_gain','PFA_top','IDR_max','dP_orf_q95'});
    all_results = [all_results; row]; %#ok<AGROW>

    if i == 1
        writetable(row, outfile);
    else
        writetable(row, outfile, 'WriteMode', 'append', 'WriteVariableNames', false);
    end
end

disp(all_results);
