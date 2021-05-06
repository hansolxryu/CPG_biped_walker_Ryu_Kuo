function [wstar,cnvrg] = gradsearch(w, x0, varargin)
% GRADSEARCH Finds a fixed point for walking model by doing
% a gradient search.
% [wstar, cnvrg] = gradsearch(w, x0) takes the model w, and returns the
% fixed point in model wstar, which includes xstar
% cnvrg indicates whether the search worked or not
% if the search fails, it returns w
% The following options can be used:
% 'delta' (1e-8) is the perturbation size for evaluating gradients
% 'info' (1) determines whether progress is printed to the screen
%            info = 1 normal, info = 2 detailed
% 'criterion' (1e-10) is the criterion for stopping the search
% 'stepsize'  (1) is the size of steps in Newton's method, which can
%             be scaled to less than one to slow the search
% 'tryharder' (1) uses multiple step sizes and also tries reversing
%             delta to find tricky fixed points

% Changes
%   v 1.1 Added 'tryharder' option, which switches on a multi-step-size
%         and reversing delta strategy for use in finding tricky fixed
%         points for the first time.
%   v 1.2 Changed a pause command to drawnow, which is the preferred way
%         to allow ctrl-c interrupts. (Art 7/24/2008)
%   v 1.3 Added provision for onestep to plot the trajectory, so make it
%         easier to track progress, requires info >= 2. Also somewhere
%         along the line somebody added tryharder.


delta = 1e-7; info = 1; criterion = 1e-9; stepsize = 1; parmvary = []; tryharder = 0;
method = 'normal';

% If walk object has a property "ddelta" for the default delta size, e.g. because of the limitations of ADAMS, use it.  Same for criterion.
try delta = get(w,'ddelta');  catch; end;  try criterion = get(w,'dcrit'); catch; end

cnvrg = 0;

if nargin < 2 || (exist('x0') && isempty(x0)) % x0 not given, get it from xstar
  if ~isempty(get(w, 'xstar'))
    x0 = get(w, 'xstar');
  else
    error('gradsearch requires an initial guess x0')
  end
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
    case 'tryharder'
      tryharder = val;
    case {'extrapolate','minparmdiff'} 
      % do nothing, option for findgait
    case {'AbsTol','RelTol', 'tf'} 
      % options for onestep
      onestepoptions = {onestepoptions{:}, opt, val};        
    otherwise
      warning('Gradsearch options: delta, info, criterion, stepsize, parmvary');
  end
end

N = get(w, 'N');
xs = x0;
wstar = w;

if info > 1
  fprintf(1,'gradsearch target: ');
  if isempty(parmvary)
    display(w);
  else
    fprintf(1, '%s = %g\n', parmvary, get(w, parmvary))
  end
  hold on;
end

enew = 0; e = 1; scount = 0; ostepsize = stepsize; calc = 1;
while max(abs(e)) > criterion
  drawnow;  % leave place for ctrl-c interrupt
  if calc
    y = zeros(N,N);
    y1 = onestep(w, xs, onestepoptions{:})-xs;
    for i=1:N
      x = xs; x(i) = x(i) + delta;
      y(:,i) = (onestep(w, x, onestepoptions{:})-x)'-y1';
      if info > 1
        onestep(w, x, onestepoptions{:});
      end
    end
    J = y / delta;
    dx = J \(-y1');
  end
  xsn = xs + dx'*stepsize;
  enew = onestep(w, xsn, onestepoptions{:})-xsn;
  if max(abs(enew)) > max(abs(e))
    if stepsize <= ostepsize*2^(-5) | ~tryharder
      if delta < 0 | ~tryharder
        cnvrg = 0;
        if info > 0
            fprintf(1,' could not converge\n');
        end
        %xs = xsn - dx'*stepsize;
        break
      else
        % if info > 0; fprintf(1, '  e = %g\n', max(abs(e))); end
        if info > 0; fprintf(1,'  trying a negative delta\n'); end
        delta = -delta;
        stepsize = ostepsize;
        calc = 1;
      end
    else
      % if info > 0; fprintf(1, '  e = %g\n', max(abs(e))); end
      if info > 0; fprintf(1,'  trying a smaller step size\n'); end
      stepsize = stepsize/2;
      calc = 0;
    end
  else
    cnvrg = 1;
    e = enew;
    xs = xsn;
    calc = 1;
    stepsize = ostepsize;
    delta = abs(delta);
  end
  if info > 1
    fprintf(1, '  e = %g\n', max(abs(e)));
  end
  % store it
  wstar = set(wstar, 'xstar', xs); 
  scount = scount + 1;
end % while

if info > 0
  fprintf(1,'gradsearch return  %d steps cnvrg = %d\n', scount, cnvrg);
end
