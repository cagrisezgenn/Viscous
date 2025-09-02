﻿function Sa_vec = sdof_PSA_vec_aug_ode(t, ag, T_vec, zeta)
% Ã‡oklu periyotlar iÃ§in tek ODE Ã§Ã¶zÃ¼mÃ¼
% t: zaman vektÃ¶rÃ¼, ag: yer ivmesi, T_vec: periyotlar, zeta: sÃ¶nÃ¼m oranÄ±

    T_vec = T_vec(:); Np = numel(T_vec);

    % --- Ã–nce Ã¶nbelleÄŸi kontrol et ------------------------------------
    key     = psa_hash(t, ag, T_vec, zeta);
    Sa_vec  = psa_cache('get', key);
    if ~isempty(Sa_vec), return; end

    % Hesaplama
    wn = 2*pi./max(T_vec,eps);
    agf = griddedInterpolant(t, ag, 'linear','nearest');

    z0 = zeros(2*Np,1);
    dt = median(diff(t));
    opts = odeset('RelTol',1e-4,'AbsTol',1e-6,...
                  'MaxStep',max(10*dt,1e-3),'InitialStep',max(0.25*dt,1e-4));
    odef = @(tt,zz) aug_rhs(tt,zz,wn,zeta,agf);
    sol  = ode23tb(odef, [t(1) t(end)], z0, opts);

    Z = deval(sol, t).';                 % Nt x (2*Np)
    X = Z(:,1:2:end);                    % Nt x Np
    V = Z(:,2:2:end);                    % Nt x Np
    Nt = size(Z,1);
    Aabs = abs( -2*zeta*(ones(Nt,1)*wn.').*V ...
                - (ones(Nt,1)*(wn.'.^2)).*X ...
                - ag(:)*ones(1,Np) );    % Nt x Np
    Sa_vec = max(Aabs,[],1).';

    % Sonucu Ã¶nbelleÄŸe yaz
    psa_cache('set', key, Sa_vec);
end

function dz = aug_rhs(tt,zz,wn,zeta,agf)
% YardÄ±mcÄ± saÄŸ taraf fonksiyonu
    x = zz(1:2:2*numel(wn)); v = zz(2:2:2*numel(wn));
    agt = agf(tt);
    ax  = v;
    av  = -2*zeta*wn.*v - (wn.^2).*x - agt;
    dz  = zeros(2*length(wn),1);
    dz(1:2:end) = ax;
    dz(2:2:end) = av;
end

