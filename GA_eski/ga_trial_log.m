function ga_trial_log(action, payload)
%GA_TRIAL_LOG Append GA evaluation rows to a CSV file.
%  Header automatically expands according to GA set sizes returned by ga_get_bounds.%  ga_trial_log('init_header', S) -> S.file (optional) sets file path.
%  ga_trial_log('append', S) where S has fields:
%    .design_set (scalar), .ga_x [1xN] for the active set (N depends on set), .vars32 (struct with fields listed in header),
%    .zeta_eq, .okP, .okQ, .okDp, .Penalty, .J, .J1, .J2, .score,
%    .pen (struct: spring_tau, spring_slender, stroke, force_cap, dp_quant, thermal_dT, cav_frac, qsat_margin, fail_bigM),
%    .acc3_max_abs, .x10_max_abs, .acc3_area_pos, .acc3_area_neg_abs, .x10_area_pos, .x10_area_neg_abs
%  Always appends; never truncates existing logs.

    persistent FILE
    if isempty(FILE), FILE = fullfile('out','trial_log.csv'); end
    if nargin>=2 && isstruct(payload) && isfield(payload,'file'), FILE = payload.file; end
    if nargin==0, action=''; end

    switch lower(string(action))
        case 'init_header'
            if ~exist('out','dir'), try, mkdir('out'); end, end
            % If file exists but header mismatches, rotate old and write fresh header
            if exist(FILE,'file')
                try
                    fid = fopen(FILE,'r');
                    line = ''; if fid>0, line = fgetl(fid); fclose(fid); end
                catch
                    line = '';
                end
                exp = header_line();
                if ~ischar(line) || ~strcmp(strtrim(line), strtrim(exp))
                    try
                        [pth,base,ext] = fileparts(FILE);
                        ts = datestr(now,'yyyymmdd_HHMMSS');
                        backup = fullfile(pth, sprintf('%s_old_%s%s', base, ts, ext));
                        movefile(FILE, backup, 'f');
                    catch
                        % ignore rotation failure
                    end
                    write_header(FILE);
                end
            else
                write_header(FILE);
            end
        case 'append'
            S = payload; %#ok<NODEF>
            if ~exist('out','dir'), try, mkdir('out'); end, end
            row = build_row_vector(S);
            append_csv(FILE, row);
        otherwise
            % no-op
    end
end

function write_header(fp)
    h = header_cells();
    fid = fopen(fp,'w');
    fprintf(fid, '%s\n', strjoin(h, ','));
    fclose(fid);
end

function s = header_line()
    h = header_cells();
    s = strjoin(h, ',');
end

function h = header_cells()
    h = {};
      % GA design variables grouped by set with variable column counts
    set_sizes = get_set_sizes();
    for s = 1:numel(set_sizes)
        for i = 1:set_sizes(s)
            h{end+1} = sprintf('set%d_x%d', s, i); %#ok<AGROW>
        end
    end
    % Full named variables snapshot
    varnames = {'d_o','Lori','d_w','D_m','cav_sf','resFactor','Lgap','n_turn','Cd0','p_exp','Lh','CdInf','K_leak','mu_ref','Vmin_fac','hA_os','Rec','Kd','beta0','Dp','b_mu','dP_cap','F_cap','antoine_A','antoine_B','antoine_C','T0_C','Ts0_C','T_env_C','hA_o_env','hA_s_env','cp_oil','cp_steel','rho_ref','T_ref_C','alpha_rho','b_beta','pf_tau','pf_gain','pf_t_on','p_amb','mu_min_phys','softmin_eps'};
    for i=1:numel(varnames), h{end+1}=['var_' varnames{i}]; end %#ok<AGROW>
    % Metrics / penalties
    meta = {'zeta_eq','okP','okQ','okDp','Penalty','J','J1','J2','score', ...
            'pen_spring_tau','pen_spring_slender','pen_stroke','pen_force_cap', ...
            'pen_dp_quant','pen_thermal_dT','pen_cav_frac','pen_qsat_margin','pen_fail_bigM', ...
            'acc3_max_abs','x10_max_abs','acc3_area_pos','acc3_area_neg_abs','x10_area_pos','x10_area_neg_abs'};
    h = [h, meta];
end

function v = build_row_vector(S)
    pens = [getfield_default(S.pen,'spring_tau',0), ...
            getfield_default(S.pen,'spring_slender',0), ...
            getfield_default(S.pen,'stroke',0), ...
            getfield_default(S.pen,'force_cap',0), ...
            getfield_default(S.pen,'dp_quant',0), ...
            getfield_default(S.pen,'thermal_dT',0), ...
            getfield_default(S.pen,'cav_frac',0), ...
            getfield_default(S.pen,'qsat_margin',0), ...
            getfield_default(S.pen,'fail_bigM',0)];
    tail = [S.zeta_eq, double(S.okP), double(S.okQ), double(S.okDp), S.Penalty, S.J, S.J1, S.J2, S.score, pens, ...
            S.acc3_max_abs, S.x10_max_abs, S.acc3_area_pos, S.acc3_area_neg_abs, S.x10_area_pos, S.x10_area_neg_abs];
     % Design variables snapshot distributed into their GA sets
     % Design variables snapshot distributed into their GA sets
    set_sizes = get_set_sizes();
    total = sum(set_sizes);
    wset = nan(1,total);
    % Try to populate design variables for all sets using vars32 snapshot
    if isfield(S,'vars32') && isstruct(S.vars32)
        set_names = {};
        for s = 1:numel(set_sizes)
            [~,~,~,nm] = ga_get_bounds(s);
            set_names = [set_names, nm]; %#ok<AGROW>
        end
        for i = 1:min(numel(set_names), total)
            nm = set_names{i};
            if isfield(S.vars32, nm)
                wset(i) = S.vars32.(nm);
            end
        end
    end
    % Overwrite the active set with ga_x values if provided
    ds = getfield_default(S,'design_set',nan);
    if isfinite(ds) && ds>=1 && ds<=numel(set_sizes)
        gx = getfield_default(S,'ga_x',[]);
        n = min(numel(gx), set_sizes(ds));
        startIdx = 1 + sum(set_sizes(1:ds-1));
        wset(startIdx:startIdx+n-1) = gx(1:n);
    end
    % vars32 in fixed order
    order = {'d_o','Lori','d_w','D_m','cav_sf','resFactor','Lgap','n_turn','Cd0','p_exp','Lh','CdInf','K_leak','mu_ref','Vmin_fac','hA_os','Rec','Kd','beta0','Dp','b_mu','dP_cap','F_cap','antoine_A','antoine_B','antoine_C','T0_C','Ts0_C','T_env_C','hA_o_env','hA_s_env','cp_oil','cp_steel','rho_ref','T_ref_C','alpha_rho','b_beta','pf_tau','pf_gain','pf_t_on','p_amb','mu_min_phys','softmin_eps'};
    w32 = nan(1,numel(order));
    if isfield(S,'vars32') && isstruct(S.vars32)
        for i=1:numel(order)
            if isfield(S.vars32, order{i})
                w32(i) = S.vars32.(order{i});
            end
        end
    end
    v = [wset, w32, tail];
    v(~isfinite(v)) = 0;
end

function append_csv(fp, row)
    fid = fopen(fp,'a');
    fmt = [repmat('%.16g,',1,numel(row)-1), '%.16g\n'];
    fprintf(fid, fmt, row);
    fclose(fid);
end



function sizes = get_set_sizes()
    persistent SIZES
    if isempty(SIZES)
        SIZES = zeros(1,5);
        for s = 1:5
            [lb,~,~,~] = ga_get_bounds(s);
            SIZES(s) = numel(lb);
        end
    end
    sizes = SIZES;
end