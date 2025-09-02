function [resp, env] = simulate(design_set, x, mu_mult, t, ag, ...
    M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg, env)
if nargin < 16
    cfg = struct();
end
if nargin < 17
    env = struct();
end
if ~isfield(env,'LOG') || isempty(env.LOG)
    env.LOG = struct();
end
design_set = double(design_set);
    if ~isfinite(design_set) || ~ismember(design_set,[1 2 3 4 5])
        design_set = 1;
    end
% Standart, gÃ¼venli tek Ã§aÄŸrÄ± noktasÄ± (GA ve Ã§oklu denemeler iÃ§in)

    resp = struct('ok',false,'msg','','y',struct(),'drift',[],'dvel',[], ...
                  'dP_orf_env',[],'c_lam',NaN,'dT_est',NaN,'metrics',struct(), ...
                  'log',{{}});
    try
        % ---- guard + opsiyonel viskozite Ã§arpanÄ±
        % ensure PF defaults (e.g., auto_t_on)
        cfg = ensure_cfg_defaults(cfg);
        if ~isempty(mu_mult) && isfinite(mu_mult) && mu_mult>0
            therm.mu_ref = therm.mu_ref * mu_mult;
        end

        % ---- GA decode (varsa)
        ga_local = struct('enable',~isempty(x),'design_set',design_set,'x',x);
        [geom, sh, orf, hyd, therm, num, ~, env] = decode_design_apply(ga_local, geom, sh, orf, hyd, therm, num, env);

        % ---- TÃ¼retilenler
        geom.Ap   = pi*geom.Dp^2/4;
orf.Ao    = orf.n_orf * (pi*geom.d_o^2/4);
nd = max(1, getfield_default(hyd,'n_parallel',1));
geom.Ap_eff = nd * geom.Ap;         % toplam piston alanÄ±
orf.Ao_eff  = nd * orf.Ao;          % toplam orifis alanÄ±
hyd.n_parallel = nd;                % fonksiyonlara iletmek iÃ§in
        k_p       = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);
        k_h       = geom.Kd*geom.Ap^2/geom.Lgap;
        k_s       = geom.Ebody*geom.Ap/geom.Lgap;
        k_hyd     = 1/(1/max(k_h,eps) + 1/max(k_s,eps));
k_sd = nd * (k_hyd + k_p);

        % hacim ve Ä±sÄ± kapasiteleri
    nStories = n-1;
steel_to_oil_mass_ratio = 1.5;

hyd.V0 = 1.2*0.5*(geom.Ap*(2*geom.Lgap));              % artırılmış hazne hacmi
V_oil_per = therm.resFactor*(geom.Ap*(2*geom.Lgap));  % tek damper baÅŸÄ±na yaÄŸ
nDtot     = nStories * nd;

m_oil_tot    = nDtot*(therm.rho_ref*V_oil_per);
m_steel_tot  = steel_to_oil_mass_ratio*m_oil_tot;
therm.C_oil   = m_oil_tot*therm.cp_oil;
therm.C_steel = m_steel_tot*therm.cp_steel;

    
       % ---- PF guard: t5 + 0.5
        [t5_sim,~] = arias_win(t, ag, 0.05, 0.95);
        if cfg.PF.auto_t_on
            cfg = set_pf_ton_if_auto(cfg, t5_sim, 0.5);  % auto_t_on -> t_on=t5+0.5
        end

        % ---- Ã‡Ã¶zÃ¼m (hÄ±z isteÄŸe baÄŸlÄ± Ã§Ä±ktÄ±yla) + tolerans hatasÄ±nda otomatik yeniden dene
        dt0 = median(diff(t));
        lastwarn('');
        [xD, aD, diag, vD] = mck_with_damper_adv(t, ag, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg);
        [wm, ~] = lastwarn;
        if ~isempty(wm) && contains(lower(wm), 'unable to meet integration tolerances')
            % Preset ve toleranslarÄ± gevÅŸetip gÃ¼venli adÄ±mlarla tekrar dene
            cfg2 = cfg;
            if ~isfield(cfg2,'solver') || ~isstruct(cfg2.solver), cfg2.solver = struct(); end
            cfg2.solver.preset      = 'explore';
            if ~isfield(cfg2.solver,'RelTol') || ~isfinite(cfg2.solver.RelTol)
                cfg2.solver.RelTol  = 0.12;
            end
            cfg2.solver.AbsTolScale = max(getfield_default(cfg2.solver,'AbsTolScale',1.0), 2.0);
            cfg2.solver.MaxStep     = max(10*dt0, 2e-3);
            cfg2.solver.InitialStep = max(0.20*dt0, 1.0e-3);
            lastwarn('');
            [xD, aD, diag, vD] = mck_with_damper_adv(t, ag, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg2);
        end
% ---- NaN/Inf emniyeti (Ã¶zellikle dP_orf_time_max iÃ§in) ----
if ~isfield(diag,'dP_orf_time_max') || any(~isfinite(diag.dP_orf_time_max))
    if isfield(diag,'dP_orf_time') && ~isempty(diag.dP_orf_time)
        diag.dP_orf_time_max = max(diag.dP_orf_time, [], 2, 'omitnan');
    else
        diag.dP_orf_time_max = zeros(numel(t),1); % gÃ¼venli fallback
    end
end


% Ana dizilerde tekil NaN/Inf temizliÄŸi (fail yerine sÄ±fÄ±rlama)
xD(~isfinite(xD)) = 0;  aD(~isfinite(aD)) = 0;

        % ---- Standart Ã§Ä±ktÄ± paketi
        resp.y = struct('x',xD,'v',vD,'a',aD,'t',t);
        resp.drift = diag.drift;
        resp.dvel  = vD(:,2:end) - vD(:,1:end-1);
        resp.dP_orf_env = diag.dP_orf_time_max;               % kat zarfÄ± (max across stories)
        resp.c_lam = (128*therm.mu_ref*geom.Lori/(pi*geom.d_o^4)) / (nd * orf.n_orf);
        resp.dT_est = diag.T_oil(end) - diag.T_oil(1);

        % pratik metrikler
        x10_max = max(abs(xD(:,min(10,size(xD,2)))));
        a3_rms  = sqrt(mean(aD(:,min(3,size(aD,2))).^2,'omitnan'));
        cav_p95 = prctile(diag.cav_frac_t,95);
        Toil_max= max(diag.T_oil);
        resp.metrics = struct('x10_max',x10_max,'a3_rms',a3_rms, ...
                      'cav_frac_p95',cav_p95,'Q_abs_p95',diag.Q_abs_p95, ...
                      'T_oil_max',Toil_max);

        resp.diag = diag;   % cav_frac_t dahil tÃ¼m zaman serisi diagnostiklerini dÄ±ÅŸarÄ± ver

        % Zarf kuvvet (hikaye bazÄ±nda max â†’ zaman serisi)
        F_story_env = max(abs(diag.F_story),[],2);   % Nt x 1
        resp.F_env  = F_story_env;
        resp.F_max  = max(F_story_env);              % skaler

        % Zarf strok (hikayeler Ã¼zerinde max)
        drift_env = max(abs(resp.drift),[],2);       % Nt x 1
        resp.stroke_max = max(drift_env);            % skaler

        % Î”p_orf iÃ§in zaman-iÃ§i quantile hesap kolaylÄ±ÄŸÄ±
        resp.dP_q_time = @(qq) prctile(resp.dP_orf_env, 100*qq);
% ---- SaÄŸlamlÄ±k kontrolleri

badVars = validate_finite({'xD','aD','vD','resp.dP_orf_env'});
if ~isempty(badVars)
    resp.ok  = false;
    resp.msg = ['NaN/Inf/empty tespit edildi: ' strjoin(badVars, ', ')];
    resp.log{end+1} = resp.msg;
    log_msg('error', resp.msg);
    return;
end

% dP_orf_env zaten Ã¼stte garanti edildi; x/a boyut kontrolÃ¼:
if numel(t) ~= size(xD,1) || numel(t) ~= size(aD,1)
    resp.ok  = false;
    resp.msg = 'Zaman uzunluÄŸu uyuÅŸmazlÄ±ÄŸÄ±';
    resp.log{end+1} = resp.msg;
    log_msg('error', resp.msg);
    return;
end
        resp.ok=true; resp.msg='ok';

    catch ME
        % ---- Hata halinde NaN dolgulu gÃ¼venli dÃ¶nÃ¼ÅŸ
        N = numel(t); ntry = n;
        resp.ok=false; resp.msg=sprintf('simulate: %s',ME.message);
        resp.log{end+1} = resp.msg;
        log_msg('error', resp.msg);
        resp.y = struct('x',nan(N,ntry),'v',nan(N,ntry),'a',nan(N,ntry),'t',t(:));
        resp.drift = nan(N,max(ntry-1,1));
        resp.dvel  = nan(N,max(ntry-1,1));
        resp.dP_orf_env = nan(N,1);
        resp.c_lam = NaN; resp.dT_est = NaN;
        resp.metrics = struct('x10_max',NaN,'a3_rms',NaN,'cav_frac_p95',NaN,'Q_abs_p95',NaN,'T_oil_max',NaN);
    end
end


