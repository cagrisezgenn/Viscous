function [updates, T] = sizing_param_diff(P_old, P_sized, gainsPF, sizing)
% Raporda: sadece DEĞİŞEN parametreler (old -> new)

% --- güvenli getter ---
function v = getf(s,f,def)
    if nargin<3, def = []; end
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = def;
    end
end

% --- Eski değerler (üst seviye veya orf.*) ---
old_d_o   = getf(P_old,'d_o',  getf(getf(P_old,'orf',struct()),'d_o',NaN));
old_Lori  = getf(P_old,'Lori', getf(getf(P_old,'orf',struct()),'L_orif',NaN));
old_n_orf = getf(P_old,'n_orf',NaN);
if isfinite(old_n_orf) && isfinite(old_d_o)
    old_Ao = old_n_orf*pi*(old_d_o^2)/4;
else
    old_Ao = getf(P_old,'A_o',NaN);
end
old_Qcap  = getf(P_old,'Qcap_big',NaN);
old_tau   = getf(getf(getf(P_old,'cfg',struct()),'PF',struct()),'tau',NaN);
old_gain  = getf(getf(getf(P_old,'cfg',struct()),'PF',struct()),'gain',NaN);
old_tg    = getf(P_old,'toggle_gain',[]);

% --- Yeni/önerilen değerler ---
new_d_o   = getf(getf(P_sized,'orf',struct()),'d_o', getf(P_sized,'d_o',NaN));
new_Lori  = getf(getf(P_sized,'orf',struct()),'L_orif', getf(P_sized,'Lori', getf(sizing,'L_orif',NaN)));
new_n_orf = getf(P_sized,'n_orf',NaN);
if isfinite(new_n_orf) && isfinite(new_d_o)
    new_Ao = new_n_orf*pi*(new_d_o^2)/4;
else
    new_Ao = getf(P_sized,'A_o',NaN);
end
new_Qcap  = getf(P_sized,'Qcap_big',NaN);
new_tau   = getf(gainsPF,'PF_tau',NaN);
new_gain  = getf(gainsPF,'PF_gain',NaN);

% --- toggle_gain vektörü (alt=3 kat, üst=2 kat) ---
try
    nStories = size(P_old.M,1)-1;
catch
    nStories = numel(old_tg); if isempty(nStories) || nStories==0, nStories = 10; end
end
tg = ones(nStories,1)*getf(gainsPF,'g_mid',1);
loN = min(3,nStories); if loN>0, tg(1:loN) = getf(gainsPF,'g_lo',1); end
hiN = min(2,nStories); if hiN>0, tg(end-hiN+1:end) = getf(gainsPF,'g_hi',1); end
new_tg = tg;

% --- skaler eşitlik denetleyici (taşma güvenli) ---
function tf = eqnum(a,b)
    if ~(isnumeric(a) && isnumeric(b)) || isempty(a) || isempty(b)
        tf = false; return;
    end
    a1 = a(1); b1 = b(1); % skalerleştir
    if ~(isfinite(a1) && isfinite(b1))
        tf = false; return;
    end
    tol = max(1e-12, 1e-6*max(1, max(abs([a1 b1]))));
    tf  = abs(a1-b1) <= tol;
end

changed = @(oldv,newv) ~( (isnumeric(oldv)&&isnumeric(newv)&&eqnum(oldv,newv)) || isequal(oldv,newv) );

rows = {
 'd_o (m)',           old_d_o,   new_d_o,   'd_o = %.6g;';
 'Lori (m)',          old_Lori,  new_Lori,  'Lori = %.6g;';
 'n_orf (-)',         old_n_orf, new_n_orf, 'n_orf = %d;';
 'A_o (m^2)',         old_Ao,    new_Ao,    'A_o = %.6g;';
 'Qcap_big (m^3/s)',  old_Qcap,  new_Qcap,  'Qcap_big = %.6g;';
 'cfg.PF.tau (s)',    old_tau,   new_tau,   'cfg.PF.tau = %.6g;';
 'cfg.PF.gain (-)',   old_gain,  new_gain,  'cfg.PF.gain = %.6g;'
};

Name = {}; Old = []; New = []; Template = {};
for i=1:size(rows,1)
    if changed(rows{i,2}, rows{i,3})
        Name{end+1,1} = rows{i,1}; %#ok<AGROW>
        Old(end+1,1)  = rows{i,2}; %#ok<AGROW>
        New(end+1,1)  = rows{i,3}; %#ok<AGROW>
        Template{end+1,1} = rows{i,4}; %#ok<AGROW>
    end
end
T = table(Name, Old, New, Template);

% --- toggle_gain farkı ---
tg_changed = true;
if ~isempty(old_tg) && numel(old_tg)==numel(new_tg)
    tg_changed = any(abs(old_tg(:)-new_tg(:))>1e-6);
end
if tg_changed
    T = [T; table({"toggle_gain ((n-1)x1)"}, NaN, NaN, {''}, 'VariableNames', T.Properties.VariableNames)];
end



% --- CSV ---
try
    outdir = fullfile('out'); if ~exist(outdir,'dir'), mkdir(outdir); end
    writetable(T, fullfile(outdir,'sizing_updates.csv'));
catch
end

% --- parametre.m’e yapıştırmalık ---
updates = struct(); updates.lines = {};
for i=1:height(T)
    if T.Template{i}~=""
        if contains(T.Name{i},"n_orf")
            updates.lines{end+1} = sprintf(T.Template{i}, round(T.New(i))); %#ok<AGROW>
        else
            updates.lines{end+1} = sprintf(T.Template{i}, T.New(i)); %#ok<AGROW>
        end
    end
end
updates.toggle_gain_vec = new_tg;
end
