function [wstar,cnvrg] = gradsearchspeedstep(w, x0, speed0, steplength0, varargin)
% GRADSEARCHSPEEDSTEP Finds a fixed point for walking model by doing
% a gradient search to match a desired speed and step length.
% [wstar, cnvrg] = gradsearchspeedstep(w, x0, speed0, steplen0) takes the model w, and returns the
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
% see also gradsearchspeed, which only tries to match speed by varying
% a single parameter.

% Changes
%   v 1.1 Added 'tryharder' option, which switches on a multi-step-size
%         and reversing delta strategy for use in finding tricky fixed
%         points for the first time.
%   v 1.2 Changed a pause command to drawnow, which is the preferred way
%         to allow ctrl-c interrupts. (Art 7/24/2008)
%   v 2.0 Made this into gradsearchspeed (Art 10/2010). Note that this may not
%         works on all models, because gaitspeed may need to be modified.
%         See walksw2 for example of mods.

delta = 1e-7; info = 1; criterion = 1e-8; stepsize = 1; parmvary1 = 'gamma'; tryharder = 1;
parmvary2 = 'Kp';,
method = 'normal';

cnvrg = 0; w0 = w; % save original w

if nargin < 2 || (exist('x0','var') && isempty(x0)) % x0 not given, get it from xstar
  if ~isempty(get(w, 'xstar'))
    x0 = get(w, 'xstar');
  else
    error('gradsearchspeedstep requires an initial guess x0')
  end
end

if nargin < 3 || isempty(speed0) % supply a speed
  speed0 = 0.4;
end

if nargin < 4 || isempty(steplength0) % supply a step length
  steplength0 = speed0 / (1.8/sqrt(9.81));
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
    case 'parmvary1'
      parmvary1 = val;
    case 'parmvary2'
      parmvary2 = val;
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
      warning('Gradsearchspeedstep options: delta, info, criterion, stepsize, parmvary');
  end
end

N = length(x0);  % base N on the number of states supplied
xs = x0; 
wstar = w;

if info >= 1
  fprintf(1,'gradsearchspeedstep target: ');
  if isempty(parmvary1)
    display(w);
  else
    fprintf(1, '%s = %g, %s = %g\n', parmvary1, get(w, parmvary1),...
      parmvary2, get(w,parmvary2))
  end
end

parms = get(w, 'parms'); 
theparameter1 = getfield(parms, parmvary1);
theparameter2 = getfield(parms, parmvary2);
%gamma = parms.gamma; 

enew = 0; e = 1; scount = 0; ostepsize = stepsize; calc = 1;
while max(abs(e)) > criterion
  drawnow;  % leave place for ctrl-c interrupt
  if calc
    y = zeros(N+2,N+2);
    y1 = onestepspeedstep(w, xs);
    for i=1:N
      x = xs; x(i) = x(i) + delta;
      y(:,i) = onestepspeedstep(w, x)'-y1';
    end
    %xc2 = onestepspeed(set(w, 'gamma', delta),xs);
    xc2 = onestepspeedstep(set(w, parmvary1, theparameter1+delta),xs);
    y(:,N+1) = xc2'  - y1'; 
    xc3 = onestepspeedstep(set(w, parmvary2, theparameter2+delta),xs);
    y(:,N+2) = xc3'  - y1'; 
    J = y / delta;
    dx = J \(-y1');
  end
  dparm1 = dx(N+1)*stepsize;
  dparm2 = dx(N+2)*stepsize;
  xsn(1:N) = xs(1:N) + dx(1:N)'*stepsize; 
  %wtry = set(w,'gamma', gamma+dgamma);
  wtry = set(w,parmvary1, theparameter1 + dparm1, parmvary2, theparameter2+dparm2);
  enew = onestepspeedstep(wtry, xsn);
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
    w = wtry;
    theparameter1 = theparameter1 + dparm1; %gamma = gamma+dgamma;
    theparameter2 = theparameter2 + dparm2;
    calc = 1;
    stepsize = ostepsize;
    delta = abs(delta);
  end
  if info > 1
    fprintf(1, '  e = %g\n', max(abs(e)));
  end
  % store results in wstar
  wstar = set(wstar, 'xstar', xs, parmvary1, theparameter1, parmvary2, theparameter2); 
  scount = scount + 1;
end

if info > 0
  fprintf(1,'gradsearch return  %d steps cnvrg = %d\n', scount, cnvrg);
end

% here is a nested function, returns error from desired gait, including
% speed as an output
function xout = onestepspeedstep(wnew, xsnew)
  
  [xc,tc] = onestep(wnew, xsnew, onestepoptions{:});
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
  xout = [xc-xsnew speed-speed0 steplength-steplength0];
  %end
end

end % main function
