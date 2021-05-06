function [wstar,cnvrg,sol] = bvpsearch(walk, x0, varargin)
% BVPSEARCH Finds a fixed point for walking model by doing
% a boundary value problem collocation search.
% [wstar, cnvrg,sol] = bvpsearch(w, x0) takes the model walk, and returns the
% fixed point in model wstar, which includes xstar
% cnvrg indicates whether the search worked or not
% if the search fails, it returns w
% The output sol contains fields:
%   x - same as time, except normalized 0 to 1
%   y - the states x, arranged with time in columns (tranpose of normal)
%   parameters - the actual step period and other bvp parameters
% The following options can be used:
% 'delta' (1e-8) is the perturbation size for evaluating gradients
% 'info' (1) determines whether progress is printed to the screen
%            info = 1 normal, info = 2 detailed
% 'criterion' (1e-8) is the criterion for stopping the search
% 'stepsize'  (1) is the size of steps in Newton's method, which can
%             be scaled to less than one to slow the search
% 'NMax'      (250) Max # of mesh points to use in bvp
% 'solver' (@bvp5c2) Specify the handle of the bvp solver to use. Pre-
%          Matlab 2008b, use bvp4c instead of bvp5c

% Changes
% Added by Art 6/2009 as a way to improve upon gradsearch

% This program needs the following functions defined for a walk class:
%   bvpsetup - Sets up initial guess for solver, a mesh solution
%   bvpfwalk - Returns state-derivative for an assumed step period
%   bvpbcfun - Residual of boundary conditions, to be reduced by solver

delta = 1e-8; info = 1; criterion = 1e-8; stepsize = 1; parmvary = []; 
reltol = 1e-4; NMax = 250; bvpsolver = @bvp5c;

cnvrg = 0; caught = 0;

if nargin < 2 || (exist('x0') && isempty(x0)) % x0 not given, get it from xstar
  if ~isempty(get(walk, 'xstar'))
    x0 = get(walk, 'xstar');
  else
    error('bvpsearch requires an initial guess x0')
  end
end

% see if x0 is a guess at a fixed point, or whether it is a structure
% containing the guess at an entire trajectory (fields x, y, parameters)
if isnumeric(x0)
  solinit = bvpsetup(walk, x0); % must be provided for the walk class
elseif isstruct(x0) % looks like it is a structure suitable for solinit
  solinit = x0;
  x0 = solinit.y(:,1);
end
  
onestepoptions = {};
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'delta'
      delta = val;
    case 'info'
      info = val;
    case 'criterion'
      criterion = val;
    case 'stepsize'
      stepsize = val;
    case 'parmvary'
      parmvary = val;
    case 'method'
      method = val;
    case 'NMax'
      NMax = val;
    case 'solver'
      bvpsolver = val;
    case {'extrapolate','minparmdiff'} 
      % do nothing, option for findgait
    case {'AbsTol','RelTol', 'tf'} 
      % options for onestep
      onestepoptions = {onestepoptions{:}, opt, val};        
    otherwise
      warning('Bvpsearch options: delta, info, criterion, stepsize, parmvary');
  end
end

N = get(walk, 'N');
wstar = walk;

statoption = 'off';
if info > 1
  fprintf(1,'bvpsearch target: ');
  if isempty(parmvary)
    display(walk);
  else
    fprintf(1, '%s = %g\n', parmvary, get(walk, parmvary))
  end
  statoption = 'on';
end

walkparms = get(walk, 'parms');

options = bvpset('Vectorized','off','AbsTol', criterion, 'RelTol', reltol, 'Stats', statoption, 'NMax', NMax);

try
  sol = bvpsolver(@bvpfwalkparm, @bvpbcfunparm, solinit, options);
catch ME
  if isempty(findstr('SingJac', ME.identifier)), % a real error
    rethrow(ME); % anything but the singular Jacobian is just an error
  end
  % consider this the same as a failed search
  caught = 1;
  sol = solinit; % return the initial guess, unmodified
  if info >= 1, fprintf('bvp solver could not converge in bvpsearch\n'); end
end
  
% evaluate whether we met the residual
bcres = bvpbcfun(walk, sol);
cnvrg = all(abs(bcres) < criterion) && ~caught;

wstar = set(walk, 'xstar', sol.y(:,1)');

if info > 0
  fprintf(1,'bvpsearch return    cnvrg = %d  bcres = %g\n',  cnvrg, max(abs(bcres)));
end

% the following pass the parameters and walk class to the fwalk and bc
% functions
  function xdot = bvpfwalkparm(varargin)
    xdot = bvpfwalk(varargin{:}, walk, walkparms);
  end
  function bcres = bvpbcfunparm(varargin)
    bcres = bvpbcfun(varargin{:}, walk, walkparms);
  end

end % bvpsearch function
