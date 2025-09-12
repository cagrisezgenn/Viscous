﻿function Sab = sdof_PSA_band_avg_aug(t, ag, T1, zeta, band_fac, Np)
% Band ortalamasÄ± alÄ±nmÄ±ÅŸ spektral ivme
    Tvec = linspace(band_fac(1)*T1, band_fac(2)*T1, Np);

    % Ã–nce Ã¶nbellekte var mÄ± kontrol et
    key = psa_hash(t, ag, Tvec, zeta);
    Sa  = psa_cache('get', key);
    if isempty(Sa)
        Sa = sdof_PSA_vec_aug_ode(t, ag, Tvec, zeta);
    end
    Sab = mean(Sa);
end

