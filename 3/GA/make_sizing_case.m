function [sizing, P_sized, S_worst] = make_sizing_case(scaled, params, gainsPF, opts)
%MAKE_SIZING_CASE Deterministic sizing pass after GA gains/PF are fixed.
%   [SIZING, P_SIZED, S_WORST] = MAKE_SIZING_CASE(SCALED, PARAMS, GAINSPF, OPTS)
%   freezes {g_lo,g_mid,g_hi, PF_tau, PF_gain} and derives
%   {d_o, n_orf, L_orif, Cd, A_o, Qcap_big} using simple rules consistent
%   with the ODE model and constraints used during GA.
%
% Inputs
%   scaled  - ground motion set used in GA
%   params  - struct with structural + damper base parameters
%   gainsPF - struct with fields: g_lo, g_mid, g_hi, PF_tau, PF_gain
%   opts    - optional settings:
%               dp_allow_frac (default 0.60)
%               alpha_lam     (default 0.15)
%               n_orf_set     (default [5 6])
%               dmm_step      (default 0.05 mm)
%               Cd_init       (default params.orf.CdInf or 0.70)
%               mu_nom        (default params.mu_ref)
%
% Outputs
%   sizing  - derived quantities (Q95_worst, Ao, d_o, n_orf, Cd, Lori, ...)
%   P_sized - params with {n_orf, orf.d_o, A_o, Qcap_big} updated and
%             gains/PF frozen according to gainsPF
%   S_worst - output of run_batch_windowed(scaled, P_sized, O) for reporting

    % 0-arg/low-arg convenience: try to pull from workspace and latest GA CSVs
    if nargin < 3 || isempty(scaled) || isempty(params) || isempty(gainsPF)
        % scaled / params from base workspace
        try
            if nargin < 1 || isempty(scaled), scaled = evalin('base','scaled'); end
        catch
        end
        try
            if nargin < 2 || isempty(params), params = evalin('base','params'); end
        catch
        end
        % Best gains/PF from latest ga_knee.csv or ga_front.csv
        if nargin < 3 || isempty(gainsPF)
            try
                dd = dir(fullfile('out','ga_*'));
                assert(~isempty(dd), 'make_sizing_case:auto: no GA output under out/ga_*');
                [~,ix] = max([dd.datenum]);
                latestDir = fullfile(dd(ix).folder, dd(ix).name);
                kneeCsv = fullfile(latestDir,'ga_knee.csv');
                if exist(kneeCsv,'file')
                    Tbest = readtable(kneeCsv);
                    row = 2; if height(Tbest) < 2, row = 1; end  % baseline at 1, knee at 2
                else
                    frontCsv = fullfile(latestDir,'ga_front.csv');
                    assert(exist(frontCsv,'file')==2, 'make_sizing_case:auto: ga_front.csv not found');
                    Tbest = readtable(frontCsv);
                    % knee selection on-the-fly
                    f1v = Tbest.f1; f2v = Tbest.f2;
                    f1n = (f1v - min(f1v)) ./ max(eps, (max(f1v)-min(f1v)));
                    f2n = (f2v - min(f2v)) ./ max(eps, (max(f2v)-min(f2v)));
                    [~,row] = min(hypot(f1n, f2n));
                end
                req = {'g_lo','g_mid','g_hi','PF_tau','PF_gain'};
                assert(all(ismember(req, Tbest.Properties.VariableNames)), ...
                    'make_sizing_case:auto: required columns missing in GA CSV');
                gainsPF = struct('g_lo',Tbest.g_lo(row), 'g_mid',Tbest.g_mid(row), ...
                                 'g_hi',Tbest.g_hi(row), 'PF_tau',Tbest.PF_tau(row), ...
                                 'PF_gain',Tbest.PF_gain(row));
            catch ME
                error('make_sizing_case:auto_init','Failed to auto-infer gains/PF from GA outputs: %s', ME.message);
            end
        end
    end

    if nargin < 4 || isempty(opts), opts = struct(); end

    % --- handy getter with default ---
    function v = getopt(s, f, def)
        if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
            v = s.(f);
        else
            v = def;
        end
    end

    % Defaults
    dp_allow_frac = getopt(opts,'dp_allow_frac',0.60);   % Δp_kv,target = 0.6*Δp_cap
    alpha_lam     = getopt(opts,'alpha_lam',0.15);       % laminer / türbülans oran hedefi
    n_orf_set     = getopt(opts,'n_orf_set',[5 6]);      % izinli delik sayıları
    dmm_step      = getopt(opts,'dmm_step',0.05);        % d_o kuantizasyon (mm)
    Cd_init       = getopt(opts,'Cd_init',getopt(getopt(params,'orf',struct()),'CdInf',0.70));
    mu_nom        = getopt(opts,'mu_nom',getopt(params,'mu_ref',0.9));  % Pa·s (nominal)

    % 1) PF ve geometry gains’i sabitle
    P = params;
    % gains → toggle_gain (alt 3, orta, üst 2)
    nStories = size(P.M,1)-1;
    tg = ones(nStories,1) * gainsPF.g_mid;
    loN = min(3,nStories); if loN>0, tg(1:loN) = gainsPF.g_lo; end
    hiN = min(2,nStories); if hiN>0, tg(end-hiN+1:end) = gainsPF.g_hi; end
    P.toggle_gain = tg;
    % PF
    if isfield(P,'cfg') && isfield(P.cfg,'PF')
        P.cfg.PF.tau  = gainsPF.PF_tau;
        P.cfg.PF.gain = gainsPF.PF_gain;
    end

    % 2) Fixed-gain değerlendirme (IO kapalı)
    O = struct('do_export',false,'quiet',true,'thermal_reset','each','order','natural', ...
               'use_orifice',true,'use_thermal',true, ...
               'mu_factors',[0.75 1.00 1.25], 'mu_weights',[0.2 0.6 0.2], 'thr', []);
    S_worst = run_batch_windowed(scaled, P, O);

    % 3) Q95 (worst) ve Δp_cap → Δp_kv,target
    % Pull from summary table (aggregate across records)
    Q95_worst = NaN;
    try
        vars = S_worst.table.Properties.VariableNames;
        if ismember('Q_q95_worst', vars)
            v = S_worst.table.Q_q95_worst;
            Q95_worst = max(v(:));
        elseif ismember('Q_q95_w', vars)
            v = S_worst.table.Q_q95_w;
            Q95_worst = max(v(:));
        end
    catch
    end
    % Model hard-limit veya thr üzerinden kap
    dp_cap_default = 1.0e9; % Pa (hard-kill ile tutarlı)
    dp_cap = dp_cap_default;
    if isfield(O,'thr') && ~isempty(O.thr) && isfield(O.thr,'dP95_max') && ~isempty(O.thr.dP95_max)
        dp_cap = O.thr.dP95_max;
    end
    dp_kv_target = dp_allow_frac * dp_cap;

    % 4) Türbülans baskın → Ao_req (Cd≈Cd_init) ve {n_orf, d_o}
    rho = P.rho;
    Ao_req = abs(Q95_worst) / max(Cd_init,eps) * sqrt(rho / max(2*dp_kv_target,eps));  % m^2
    % n_orf adayları üzerinde en yakın kuantize d_o seç
    dgrid = @(dmm) (dmm_step * round(dmm./max(dmm_step,eps))); % mm grid
    best = struct('n_orf',NaN,'d_o',NaN,'A_o',NaN,'Cd',Cd_init,'Lori',NaN);
    best_err = inf;
    for n_orf = n_orf_set(:).'
        d_o = sqrt(4*Ao_req/(pi*n_orf));       % m
        dmm = dgrid(1000*d_o); d_o_q = dmm/1000;  % kuantize m
        A_o = n_orf * (pi*d_o_q^2/4);
        err = abs(A_o - Ao_req)/max(Ao_req,eps);
        if err < best_err
            best_err = err;
            best.n_orf = n_orf; best.d_o = d_o_q; best.A_o = A_o;
        end
    end

    % 5) Laminer kol: L_orif öyle ki Δp_lam(Q95) ≤ α·Δp_kv,target
    % R_lam,tot = (128 μ L / (π d^4)) / n_orf ; Δp_lam = R_lam,tot * Q
    L_orif = alpha_lam * dp_kv_target * (pi * best.d_o^4) * best.n_orf / max(128*mu_nom*abs(Q95_worst),eps);
    best.Lori = L_orif;

    % 6) Re→Cd(Re) 1–2 iterasyon ve Ao_req güncelle (delik başına Q)
    Cd = Cd_init; Re = NaN; %#ok<NASGU>
    for it = 1:2
        Q_hole = Q95_worst / max(best.n_orf,1);
        Re = 4*rho*abs(Q_hole) / max(pi*mu_nom*best.d_o,eps); % ≈ (4 ρ Q)/(π μ d)
        Cd = P.orf.CdInf - (P.orf.CdInf - P.orf.Cd0) / (1 + (Re/max(P.orf.Rec,eps))^P.orf.p_exp);
        Ao_req = abs(Q95_worst) / max(Cd,eps) * sqrt(rho / max(2*dp_kv_target,eps));
        % d_o’yu yeniden kuantize et (n_orf sabit; gerekirse n_orf_set genişletilebilir)
        d_o = sqrt(4*Ao_req/(pi*best.n_orf));
        dmm = dgrid(1000*d_o); best.d_o = dmm/1000; best.A_o = best.n_orf*(pi*best.d_o^2/4);
    end
    best.Cd = Cd;

    % 7) Qcap_big (emniyetli)
    Qcap_big = 1.2 * abs(Q95_worst);

    % 8) P_sized üret (params üzerine yaz)
    P_sized = P;
    P_sized.n_orf    = best.n_orf;
    P_sized.orf.d_o  = best.d_o;
    P_sized.A_o      = best.A_o;
    P_sized.Qcap_big = Qcap_big;

    % 9) Sizing paketi (rapor)
    sizing = struct();
    sizing.Q95_worst   = Q95_worst;
    sizing.dp_cap      = dp_cap;
    sizing.dp_kv_target= dp_kv_target;
    sizing.n_orf       = best.n_orf;
    sizing.d_o         = best.d_o;
    sizing.A_o         = best.A_o;
    sizing.Cd          = best.Cd;
    sizing.L_orif      = best.Lori;
    sizing.Qcap_big    = Qcap_big;
    sizing.alpha_lam   = alpha_lam;
    sizing.notes       = 'Deterministic sizing pass; PF/gains fixed.';

    % --- rapor: parametre.m içinde hangi satırlar değişecek?
try
    [updates, T] = sizing_param_diff(params, P_sized, gainsPF, sizing); %#ok<NASGU>
catch
end
end
