%% Parametre ızgarası üzerinde analiz çalıştır
% Bu betik parametreler.m içindeki tarama dizilerini alır ve
% tüm kombinasyonlar için run_one_record_windowed fonksiyonunu çağırır.

% 1) Parametreleri yükle
parametreler; %#ok<*NASGU>

% Tarama yapılacak parametre dizileri
vals = {Dp_vals(:), Lgap_vals(:), d_o_vals(:), Lori_vals(:), mu_ref_vals(:), ...
        d_w_vals(:), D_m_vals(:), n_turn_vals(:), rho_vals(:), n_orf_vals(:), ...
        orf_Cd0_vals(:), orf_CdInf_vals(:), orf_Rec_vals(:), orf_p_exp_vals(:), ...
        orf_p_amb_vals(:), orf_p_cav_eff_vals(:), orf_cav_sf_vals(:), ...
        thermal_hA_vals(:), cfg_PF_tau_vals(:), cfg_PF_gain_vals(:), ...
        n_dampers_per_story(:), toggle_gain(:)};
varNames = {'Dp','Lgap','d_o','Lori','mu_ref','d_w','D_m','n_turn','rho','n_orf', ...
            'orf_Cd0','orf_CdInf','orf_Rec','orf_p_exp','orf_p_amb','orf_p_cav_eff','orf_cav_sf', ...
            'thermal_hA_W_perK','cfg_PF_tau','cfg_PF_gain','n_dampers_per_story','toggle_gain'};

% 2) Tüm kombinasyonları üret
grid = cell(1,numel(vals));
[grid{:}] = ndgrid(vals{:});
comb = zeros(numel(grid{1}), numel(vals));
for k = 1:numel(grid)
    comb(:,k) = grid{k}(:);
end

% 3) Yer hareketi kaydını yükle (hızlı olması için ilk kayıt)
[~, scaled] = load_ground_motions(T1);
rec = scaled(1);

% Baz yapılar
thermal_base = thermal;
cfg_base = cfg;

% Çıktı dosyası
outdir = 'out'; if ~exist(outdir,'dir'), mkdir(outdir); end
outfile = fullfile(outdir, 'param_grid_results.csv');
if exist(outfile,'file'), delete(outfile); end
all_results = table();

% 4) Parametre kombinasyonları üzerinde döngü
opts = struct('use_orifice', true, 'use_thermal', true);
for i = 1:size(comb,1)
    % --- Parametreleri ayıkla ---
    Dp_i    = comb(i,1);   Lgap_i = comb(i,2);   d_o_i  = comb(i,3);
    Lori_i  = comb(i,4);   mu_i   = comb(i,5);   d_w_i  = comb(i,6);
    D_m_i   = comb(i,7);   n_turn_i = comb(i,8); rho_i  = comb(i,9);
    n_orf_i = comb(i,10);  Cd0_i = comb(i,11);   CdInf_i = comb(i,12);
    Rec_i   = comb(i,13);  p_exp_i = comb(i,14); p_amb_i = comb(i,15);
    p_cav_i = comb(i,16);  cav_sf_i = comb(i,17);
    hA_i    = comb(i,18);  PF_tau_i = comb(i,19); PF_gain_i = comb(i,20);
    ndps_i  = comb(i,21);  gain_i   = comb(i,22);

    % --- Türetilmiş büyüklükler ---
    Ap_i    = pi*Dp_i^2/4;
    k_h     = Kd*Ap_i^2/Lgap_i;
    k_s     = Ebody*Ap_i/Lgap_i;
    k_hyd   = 1/(1/k_h + 1/k_s);
    k_p     = Gsh*d_w_i^4/(8*n_turn_i*D_m_i^3);
    k_sd_i  = k_hyd + k_p;
    c_lam0_i= 12*mu_i*Lori_i*Ap_i^2/d_o_i^4;
    A_o_i   = n_orf_i * (pi*d_o_i^2/4);
    orf_i = struct('Cd0',Cd0_i,'CdInf',CdInf_i,'Rec',Rec_i,'p_exp',p_exp_i, ...
                   'p_amb',p_amb_i,'p_cav_eff',p_cav_i,'cav_sf',cav_sf_i, ...
                   'd_o',d_o_i,'veps',orf.veps);
    Qcap_big_i = max(orf_i.CdInf*A_o_i,1e-9) * sqrt(2*1.0e9/rho_i);
    thermal_i = thermal_base; thermal_i.hA_W_perK = hA_i;
    cfg_i = cfg_base; cfg_i.PF.tau = PF_tau_i; cfg_i.PF.gain = PF_gain_i;
    c_lam_min_i = max(c_lam_min_abs, c_lam_min_frac*c_lam0_i);

    params = struct('M',M,'C0',C0,'K',K,'k_sd',k_sd_i,'c_lam0',c_lam0_i, ...
        'orf',orf_i,'rho',rho_i,'Ap',Ap_i,'A_o',A_o_i,'Qcap_big',Qcap_big_i,'mu_ref',mu_i, ...
        'thermal',thermal_i,'T0_C',T0_C,'T_ref_C',T_ref_C,'b_mu',b_mu, ...
        'c_lam_min',c_lam_min_i,'c_lam_cap',c_lam_cap,'Lgap',Lgap_i, ...
        'cp_oil',cp_oil,'cp_steel',cp_steel,'steel_to_oil_mass_ratio',steel_to_oil_mass_ratio, ...
        'story_mask',story_mask,'resFactor',resFactor,'cfg',cfg_i,'story_height',story_height);
    params.n_dampers_per_story = ndps_i;
    params.toggle_gain = gain_i;

    % --- Simülasyon ---
    out = run_one_record_windowed(rec, [], params, opts);
    w = out.worst;

    % Koruma amaçlı yardımcılar
    if isfield(w,'x10_max_D'), x10 = w.x10_max_D; elseif isfield(w,'x10_pk_D'), x10 = w.x10_pk_D; else, x10 = NaN; end
    if isfield(w,'a10abs_max_D'), a10 = w.a10abs_max_D; elseif isfield(w,'a10abs_pk_D'), a10 = w.a10abs_pk_D; else, a10 = NaN; end
    dP95 = Utils.getfield_default(w,'dP_orf_q95',NaN);
    Qcap95 = Utils.getfield_default(w,'Qcap_ratio_q95',NaN);
    cav   = Utils.getfield_default(w,'cav_pct',NaN);
    Tend  = Utils.getfield_default(w,'T_oil_end',NaN);
    muend = Utils.getfield_default(w,'mu_end',NaN);
    PFp95 = Utils.getfield_default(w,'PF_p95',NaN);
    Qq50  = Utils.getfield_default(w,'Q_q50',NaN);
    Qq95  = Utils.getfield_default(w,'Q_q95',NaN);
    dPq50 = Utils.getfield_default(w,'dP_orf_q50',NaN);
    dPq95 = Utils.getfield_default(w,'dP_orf_q95',dP95);
    Tsteel= Utils.getfield_default(w,'T_steel_end',NaN);
    Eor   = Utils.getfield_default(w,'E_orifice_full',NaN);
    Estr  = Utils.getfield_default(w,'E_struct_full',NaN);
    Eratio= Utils.getfield_default(w,'E_ratio_full',NaN);
    Etot  = Eor + Estr;

    % Sonuç satırı
    param_vals = comb(i,:);
    metr_vals = [x10, a10, dP95, Qcap95, cav, Tend, muend, PFp95, Qq50, Qq95, dPq50, dPq95, Tend, Tsteel, Etot, Eor, Estr, Eratio];
    row = array2table([param_vals, metr_vals], 'VariableNames', [varNames, ...
        {'x10_max_damperli','a10abs_max_damperli','dP95_worst','Qcap95_worst', ...
         'cav_pct_worst','T_end_worst','mu_end_worst','PF_p95_worst','Q_q50_worst', ...
         'Q_q95_worst','dP_orf_q50_worst','dP_orf_q95_worst','T_oil_end_worst', ...
         'T_steel_end_worst','energy_tot_sum','E_orifice_sum','E_struct_sum','E_ratio'}]);
    all_results = [all_results; row]; %#ok<AGROW>

    if i == 1
        writetable(row, outfile);
    else
        writetable(row, outfile, 'WriteMode', 'append', 'WriteVariableNames', false);
    end
end

disp(all_results);

