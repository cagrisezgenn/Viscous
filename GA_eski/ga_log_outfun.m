function varargout = ga_log_outfun(varargin)
%GA_LOG_OUTFUN OutputFcn + helper for GA best-so-far logging.
%  Usage:
%    ga_log_outfun('init', cfg)  % set stage/save_path
%    best = ga_log_outfun('get')  % retrieve best so far (struct with x,f,gen,time,stage)
%    [state, options, optchanged] = ga_log_outfun(options, state, flag) % as GA OutputFcn

  persistent CFG BEST

  % Helper entry points
  if nargin>=1 && (ischar(varargin{1}) || (isstring(varargin{1}) && isscalar(varargin{1})))
      cmd = char(varargin{1});
      switch lower(cmd)
          case 'init'
              if nargin>=2 && isstruct(varargin{2})
                  CFG = varargin{2};
              else
                  CFG = struct('stage', NaN, 'save_path','');
              end
              BEST = struct();
              if nargout>0, varargout{1} = []; end
              return;
          case 'get'
              if isempty(BEST), BEST = struct(); end
              if nargout>0, varargout{1} = BEST; else, disp(BEST); end
              return;
          case 'stop'
              % Request graceful stop via sentinel file
              sf = local_stop_file(CFG);
              try
                  [d,~] = fileparts(sf); if ~isempty(d), mkdir(d); end
                  fid = fopen(sf,'w'); if fid>0, fwrite(fid,'stop'); fclose(fid); end
              catch
              end
              if nargout>0, varargout{1} = sf; end
              return;
          case 'clearstop'
              sf = local_stop_file(CFG);
              if exist(sf,'file'), try, delete(sf); catch, end, end
              if nargout>0, varargout{1} = sf; end
              return;
          case 'status'
              sf = local_stop_file(CFG);
              st = struct('stop_file',sf,'stop_requested',exist(sf,'file')==2,'stage',tern(isfield(CFG,'stage'),CFG.stage,NaN));
              if nargout>0, varargout{1} = st; else, disp(st); end
              return;
          case 'clear'
              CFG = []; BEST = struct(); if nargout>0, varargout{1} = []; end
              return;
      end
  end

  % GA OutputFcn signature: (options,state,flag)
  options = varargin{1}; %#ok<NASGU>
  state   = varargin{2};
  flag    = varargin{3}; %#ok<NASGU>
  optchanged = false;

  try
      % Compute best of current generation
      [fbest, i] = min(state.Score);
      xbest = state.Population(i,:);
      BEST = struct('x', xbest, 'f', fbest, ...
                    'gen', state.Generation, 'time', now, ...
                    'stage', tern(isfield(CFG,'stage') && ~isempty(CFG.stage), CFG.stage, NaN));
      % Save lightweight record if configured
      if ~isempty(CFG) && isfield(CFG,'save_path') && ~isempty(CFG.save_path)
          rec = BEST; %#ok<NASGU>
          try
              d = fileparts(CFG.save_path); if ~isempty(d) && ~exist(d,'dir'), mkdir(d); end
              save(CFG.save_path, 'rec');
          catch
              % ignore file I/O errors in OutputFcn
          end
      end
      % If sentinel exists, request GA to stop at this generation boundary
      sf = local_stop_file(CFG);
      if exist(sf,'file')==2
          state.StopFlag = 'user requested stop';
      end
  catch
      % ignore any issues during OutputFcn
  end

  varargout = {state, varargin{1}, optchanged};
end

function out = tern(cond, a, b)
  if cond, out = a; else, out = b; end
end

function sf = local_stop_file(C)
  try
      if isstruct(C) && isfield(C,'save_path') && ~isempty(C.save_path)
          sf = fullfile(fileparts(C.save_path), 'STOP_GA');
      else
          sf = fullfile('out','STOP_GA');
      end
  catch
      sf = fullfile('out','STOP_GA');
  end
end

