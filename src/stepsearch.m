function [wstar,cnvrg,states] = stepsearch(w, x0, varargin)
% STEPSEARCH Finds a stable fixed point for walking model by taking
% multiple steps.
% [wstar, cnvrg] = stepsearch(w, x0) takes the model w, and returns the
% fixed point in model wstar, which includes xstar
% cnvrg indicates whether the search worked or not
% states is a vector comprised of the start states for the series of
% successful steps
% if the search fails, it returns w

% 'info' (1) determines whether progress is printed to the screen
%            info = 1 normal, info = 2 detailed
% 'cnvgcriterion' (1e-12) is the criterion for stopping the search
info = 1; cnvgcriterion = 1e-9; errorcriterion = 1e-4; parmvary = [];
method = 'normal';

% If walk object has a property "dcrit" for the default criterion size, e.g. because of the limitations of ADAMS, use it
try cnvgcriterion = get(w,'dcrit'); catch; end

cnvrg = 0;

if nargin < 2 || (exist('x0') && isempty(x0)) % x0 not given, get it from xstar
  if ~isempty(get(w, 'xstar'))
    x0 = get(w, 'xstar');
  else
    error('stepsearch requires an initial guess x0')
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
    case 'cnvgcriterion'
      cnvgcriterion = val;
    case 'stepsize'
      stepsize = val;
    case 'parmvary'
      parmvary = val;
    case 'method'
      method = val;
    case 'extrapolate' % do nothing, option for findgait
    case {'AbsTol','RelTol', 'tf'} % options for onestep
      onestepoptions = {onestepoptions{:}, opt, val};
    otherwise
      warning('Stepsearch options: info, cnvgcriterion, parmvary')
  end
end


N = get(w, 'N');
xs = x0; states = x0;
wstar = w;

if info > 1
  fprintf(1,'stepsearch target: ');
  if isempty(parmvary)
    display(w);
  else
    fprintf(1, '%s = %g\n', parmvary, get(w, parmvary))
  end
end

enew = 0; e = 10; scount = 0; cnvrg = 1;
while max(abs(e)) > cnvgcriterion
  scount = scount + 1;
  xnew = onestep(w, xs, onestepoptions);
  enew = xnew - xs;
  if max(abs(enew)) >= max(abs(e)) + errorcriterion
    cnvrg = 0;
    if info > 0
      fprintf(1,' could not converge\n');
    end
    xs = x0; % give up, return original initial condition
    break
  end
  cnvrg = 1;
  xs = xnew;
  e = enew;
  if info > 1
    fprintf(1, '  e = %g\n', max(abs(e)));
  end
  % store it
  states(scount+1,:) = xs;
end

if info > 0
  fprintf(1,'stepsearch return  %d steps cnvrg = %d\n', scount, cnvrg);
end

if cnvrg
  wstar = set(wstar, 'xstar', xs);
end

