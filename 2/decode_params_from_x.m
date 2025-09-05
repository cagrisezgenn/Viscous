function P = decode_params_from_x(params0, x)
%DECODE_PARAMS_FROM_X Copy params and apply design vector x.
%   x = [d_o_mm, n_orf, g_lo, g_mid, g_hi, PF_tau, PF_gain]

    d_o_mm = x(1); n_orf = round(x(2));
    g_lo = x(3); g_mid = x(4); g_hi = x(5);
    PF_tau = x(6); PF_gain = x(7);

    P = params0;                       % copy
    P.orf.d_o = d_o_mm * 1e-3;         % mm -> m
    P.n_orf   = n_orf;
    P.A_o     = P.n_orf * (pi * P.orf.d_o^2 / 4);
    P.Qcap_big= max(P.orf.CdInf * P.A_o, 1e-9) * sqrt(2 * 1.0e9 / P.rho);

    % toggle_gain as 3-zone vector over stories: [1..3]=lo, [4..n-2]=mid, [n-1..n-?]=hi
    n  = size(P.M,1); nStories = n-1;
    tg = ones(nStories,1) * g_mid;
    loN = min(3, nStories); if loN > 0, tg(1:loN) = g_lo; end
    hiN = min(2, nStories); if hiN > 0, tg(end-hiN+1:end) = g_hi; end
    P.toggle_gain = tg;

    % Optional PF parameters under P.cfg.PF
    if isfield(P,'cfg') && isfield(P.cfg,'PF')
        P.cfg.PF.tau  = PF_tau;
        P.cfg.PF.gain = PF_gain;
    end
end

