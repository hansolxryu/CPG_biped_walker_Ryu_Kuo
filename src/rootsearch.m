function [wstar,cnvrg] = rootsearch(w, x0, varargin)
% ROOTSEARCH Finds a fixed point for walking model by doing
% a gradient search to match a desired speed.
% [wstar, cnvrg] = rootsearch(w, x0, speed0) takes the model w, and returns the
% fixed point in model wstar, which includes xstar, to walk at speed speed0
% cnvrg indicates whether the search worked or not
% if the search fails, it returns w
% The following options can be used:
% 'parmvary' ('gamma') is the parameter varied to reach the desired gait.
% 'delta' (1e-8) is the perturbation size for evaluating gradients
% 'info' (1) determines whether progress is printed to the screen
%            info = 1 normal, info = 2 detailed
% 'criterion' (1e-10) is the criterion for stopping the search
% 'stepsize'  (1) is the size of steps in Newton's method, which can
%             be scaled to less than one to slow the search
% 'tryharder' (1) uses multiple step sizes and also tries reversing
%             delta to find tricky fixed points
% 'rootfunction' (@onesteproot) is a function whose root is to be found
% 'roottarget' is an option for the root function, to indicate what value
%   to aim for

% Changes
%   v 1.1 Added 'tryharder' option, which switches on a multi-step-size
%         and reversing delta strategy for use in finding tricky fixed
%         points for the first time.
%   v 1.2 Changed a pause command to drawnow, which is the preferred way
%         to allow ctrl-c interrupts. (Art 7/24/2008)
%   v 2.0 Made this into rootsearch (Art 10/2010). Note that this may not
%         works on all models, because gaitspeed may need to be modified.
%         See walksw2 for example of mods.

delta = 1e-6; info = 1; criterion = 1e-8; stepsize = 1; tryharder = 0;
method = 'normal'; 

matchname = strcmp('rootfunction', fieldnames(w)); % see if there is a rootfunction
if any(matchname) % exists, so use it
  rootfunction = get(w, 'rootfunction');
else
  rootfunction = @onesteproot;
end

cnvrg = 0; w0 = w; % save original w

if nargin < 2 || (exist('x0') && isempty(x0)) % x0 not given, get it from xstar
  if ~isempty(get(w, 'xstar'))
    x0 = get(w, 'xstar');
  else
    error('rootsearch requires an initial guess x0')
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
    case 'method'
      method = val;
    case 'rootfunction'
      rootfunction = val;
    case 'tryharder'
      tryharder = val;
    case 'parmvary'
      parmvary = val;
      onestepoptions = {onestepoptions{:}, opt, val};
    case {'extrapolate','minparmdiff'} 
      % do nothing, option for findgait
    case {'AbsTol','RelTol', 'tf', 'roottarget','parmvary1','parmvary2'} 
      % options for onesteproot
      onestepoptions = {onestepoptions{:}, opt, val};        
    otherwise
      warning('rootsearch options: delta, info, criterion, stepsize, parmvary');
  end
end

N = length(x0);  % base N on the number of states supplied
xs = x0; 
wstar = w;

if info > 1
  fprintf(1,'rootsearch target: ');
  if isempty(parmvary)
    display(w);
  else
    fprintf(1, '%s = %g\n', parmvary, get(w, parmvary))
  end
end

parms = get(w, 'parms'); 

enew = 0; e = 1; scount = 0; ostepsize = stepsize; calc = 1;
while max(abs(e)) > criterion
  drawnow;  % leave place for ctrl-c interrupt
  if calc
    y = zeros(N,N);
    y1 = rootfunction(w, xs, onestepoptions{:});
    for i=1:N
      x = xs; x(i) = x(i) + delta;
      y(:,i) = rootfunction(w, x, onestepoptions{:})'-y1';
    end
    J = y / delta;
    dx = J \(-y1');
  end
  xsn = xs + dx'*stepsize; 
  if info > 1, 
    enew = rootfunction(w, xsn, 'plotstep', 1, onestepoptions{:});
  else
    enew = rootfunction(w, xsn, onestepoptions{:});   
  end
  if max(abs(enew)) > max(abs(e)) || any(isnan(enew)) % if error not reducing, try a smaller step size
    if stepsize <= ostepsize*2^(-5) | ~tryharder % try up to 5 times
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
  else % error is reducing, so let's keep the changes
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
  % store results in wstar
  wstar = set(wstar, 'xstar', xs);  
  scount = scount + 1;
end

if info > 0
  fprintf(1,'rootsearch return  %d steps cnvrg = %d\n', scount, cnvrg);
end

% here is a nested function, returns error from desired gait, including
% speed as an output
function xout = onestepspeed(wnew, xsnew, varargin)
  
  [xc,tc] = onestep(wnew, xsnew, varargin{:});
  %msgstr = lastwarn; % if things go well, return quickly
  %if strcmp(msgstr, 'Matrix is singular to working precision.')
  %  xout = NaN;
  %  return
  %else
  %end
  [speed, steplength, stepfreq] = gaitspeed(wnew, xsnew, xc, tc);
  if isempty(xc) % sometimes xc fails, so return NaN
    xc = xsnew*NaN;
  end
  xout = [xc-xsnew speed-speed0];
  %end
end

function xout = onesteproot(wnew, xsnew, varargin)
  
  [xc,tc] = onestep(wnew, xsnew, varargin{:});
  if isempty(xc) % sometimes xc fails, so return NaN
    xc = xsnew*NaN;
  end
  xout = xc-xsnew;
end

end % main function

