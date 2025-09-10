%% Parametre ızgarası üzerinde analiz çalıştır
% Bu betik parametreler.m içindeki tarama dizilerini alır ve
% tüm kombinasyonlar için run_one_record_windowed fonksiyonunu çağırır.

% 1) Parametreleri yükle
parametreler; %#ok<*NASGU>

% Tarama yapılacak parametre dizilerini oku
val_cells = cell(size(sweep_vars));
for ii = 1:numel(sweep_vars)
    val_cells{ii} = eval(sweep_vars{ii})(:);
end

% 2) Tüm kombinasyonları üret
[grids{1:numel(val_cells)}] = ndgrid(val_cells{:});
n_comb = numel(grids{1});

% 3) Yer hareketi kaydını yükle (hızlı olması için ilk kayıt)
[~, scaled] = load_ground_motions(T1);
rec = scaled(1);

% 4) Temel parametre yapısı (tarama dışındaki parametreler)
base_orf = orf; base_orf.d_o = orf.d_o(1);
base_params = struct('M',M,'C0',C0,'K',K,'k_sd',k_sd,'c_lam0',c_lam0(1), ...
    'orf',base_orf,'rho',rho,'Ap',Ap,'A_o',A_o(1),'Qcap_big',Qcap_big(1),'mu_ref',mu_ref, ...
    'thermal',thermal,'T0_C',T0_C,'T_ref_C',T_ref_C,'b_mu',b_mu, ...
    'c_lam_min',c_lam_min(1),'c_lam_cap',c_lam_cap,'Lgap',Lgap, ...
    'cp_oil',cp_oil,'cp_steel',cp_steel,'steel_to_oil_mass_ratio',steel_to_oil_mass_ratio, ...
    'story_mask',story_mask,'resFactor',resFactor,'cfg',cfg,'story_height',story_height,'d_o',d_o(1));

% Çıktı dosyası
outdir = 'out'; if ~exist(outdir,'dir'), mkdir(outdir); end
outfile = fullfile(outdir, 'param_grid_results.csv');
if exist(outfile,'file'), delete(outfile); end
all_results = table();

% 5) Parametre kombinasyonları üzerinde döngü
opts = struct('use_orifice', true, 'use_thermal', true);
for i = 1:n_comb
    params = base_params;
    row_vals = zeros(1, numel(sweep_vars));
    for j = 1:numel(sweep_vars)
        val = grids{j}(i);
        field = sweep_vars{j};
        params.(field) = val;
        row_vals(j) = val;
    end

    if isfield(params, 'd_o')
        d_o_val = params.d_o;
        params.A_o = n_orf * (pi*d_o_val^2/4);
        params.c_lam0 = 12*mu_ref*Lori*Ap^2 / d_o_val^4;
        params.orf.d_o = d_o_val;
        params.Qcap_big = max(params.orf.CdInf*params.A_o, 1e-9) * sqrt(2*1.0e9/rho);
        params.c_lam_min = max(c_lam_min_abs, c_lam_min_frac*params.c_lam0);
    end

    out = run_one_record_windowed(rec, [], params, opts);
    row = array2table([row_vals, out.metr.PFA_top, out.metr.IDR_max, out.metr.dP_orf_q95], ...
        'VariableNames', [sweep_vars, {'PFA_top','IDR_max','dP_orf_q95'}]);
    all_results = [all_results; row]; %#ok<AGROW>

    if i == 1
        writetable(row, outfile);
    else
        writetable(row, outfile, 'WriteMode', 'append', 'WriteVariableNames', false);
    end
end

disp(all_results);
