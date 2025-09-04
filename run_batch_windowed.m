function summary = run_batch_windowed(scaled, params, opts)
%RUN_BATCH_WINDOWED Çoklu kayıtları çalıştırır ve özetler.
%   summary = RUN_BATCH_WINDOWED(scaled, params, opts)
%
%   Dönen summary yapısı:
%     results      - kayıt bazlı çıktı (run_one_record_windowed)
%     table        - tablolaştırılmış sonuçlar
%     max_indices  - [iPFA iIDR]
%     qc           - basit kalite kontrol bilgileri

if nargin<3, opts = struct(); end

nrec = numel(scaled);
results(nrec) = struct();

for k = 1:nrec
    outk = run_one_record_windowed(scaled(k), params, opts);
    results(k) = outk;
    % termal reset politikası
    if isfield(outk.diag,'dT_est')
        switch opts.thermal_reset
            case 'carry'
                params.T0_C = params.T0_C + outk.diag.dT_est;
            case 'cooldown'
                hA = params.thermal.hA_W_perK;
                Cth = 1; % bilinmiyor
                params.T0_C = params.thermal.T_env_C + (params.T0_C + outk.diag.dT_est - params.thermal.T_env_C) * exp(-hA/Cth * opts.cooldown_s);
        end
    end
    fprintf('[%d/%d] rec="%s"  t5=%.2fs t95=%.2fs cov=%.2f  PFA=%.2f  IDR=%.4f\n', ...
        k,nrec,outk.name,outk.win.t5,outk.win.t95,outk.win.coverage,outk.metr.PFA_top,outk.metr.IDR_max);
end

PFA_vals = arrayfun(@(s) s.metr.PFA_top, results);
IDR_vals = arrayfun(@(s) s.metr.IDR_max, results);
[~,iPFA] = max(PFA_vals);
[~,iIDR] = max(IDR_vals);

summary.results = results;
summary.max_indices = struct('PFA',iPFA,'IDR',iIDR);

% coverage istatistikleri
covs = arrayfun(@(s) s.win.coverage, results);
summary.coverage_mean = mean(covs);
summary.coverage_min  = min(covs);

% tablo
names = {results.name}.'; scale = [results.scale].'; SaT1 = [results.SaT1].';
t5 = arrayfun(@(s) s.win.t5, results).';
t95= arrayfun(@(s) s.win.t95, results).';
PFA= PFA_vals.'; IDR= IDR_vals.';
dP95 = arrayfun(@(s) getfield_def(s.metr,'dP_orf_q95',NaN), results).';
Qcap95 = arrayfun(@(s) getfield_def(s.metr,'Qcap_ratio_q95',NaN), results).';
cav = arrayfun(@(s) getfield_def(s.metr,'cav_pct',NaN), results).';
T_end = arrayfun(@(s) getfield_def(s.metr,'T_oil_end',NaN), results).';
mu_end= arrayfun(@(s) getfield_def(s.metr,'mu_end',NaN), results).';
E_orf= arrayfun(@(s) getfield_def(s.metr,'E_orifice',NaN), results).';
E_str= arrayfun(@(s) getfield_def(s.metr,'E_struct',NaN), results).';
E_ratio= arrayfun(@(s) getfield_def(s.metr,'E_ratio',NaN), results).';
z1 = arrayfun(@(s) getfield_def(s.metr,'zeta1_hot',NaN), results).';
z21= arrayfun(@(s) getfield_def(s.metr,'z2_over_z1',NaN), results).';
coverage=covs.';
pass_fail = coverage>0.9;
summary.table = table(names,scale,SaT1,t5,t95,PFA,IDR,dP95,Qcap95,cav,T_end,mu_end,E_orf,E_str,E_ratio,z1,z21,coverage,pass_fail);

% QC
qc = struct();
qc.Qcap95 = all(Qcap95<0.5 | isnan(Qcap95));
qc.cav = all(cav==0 | isnan(cav));
qc.dP = all(dP95<5e7 | isnan(dP95));
qc.T_end = all(T_end<=75 | isnan(T_end));
qc.mu_end = all(mu_end>=0.5 | isnan(mu_end));
summary.qc = qc;

fprintf('Worst PFA: %s (%.2f m/s^2) | Worst IDR: %s (%.4f)\n', ...
    results(iPFA).name, results(iPFA).metr.PFA_top, ...
    results(iIDR).name, results(iIDR).metr.IDR_max);
fprintf('Coverage mean=%.2f min=%.2f\n', summary.coverage_mean, summary.coverage_min);
end

function v = getfield_def(s,f,def)
    if isfield(s,f), v = s.(f); else, v = def; end
end
