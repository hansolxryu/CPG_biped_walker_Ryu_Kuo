function [wstar,cnvrg,sol] = findgaitspeedb(w0, varargin)
% FINDGAITSPEEDB determines the impulsive-pushoff and hip spring needed
%   to produce a desired speed and steplength, using bvpsearch techniques.
%   [w, cnvrg] = findgaitspeedb(w0, 'speed', speed, 'steplength', steplength)
%   takes the existing gait w0, and tries to find P and Kp values to match
%   the desired speed and steplength, returning the new gait in w, with the
%   flag cnvrg indicating whether the search was successful
%   [w, cnvrg] = findgaitspeedb(w0, speed, steplength) takes the known
%   existing gait and tries to find P and Kp values to match the desired
%   speed and steplength, returning the new gait in w, with the flag
%   cnvrg indicating whether the search was successful.
%   findgaitspeedb(w0, speed, steplength, option, val, ...) applies the
%   desired options, which include the following:
%   'dPrel', 'dKprel' (0.001) relative perturbation sizes for P and Kp to
%      determine gradients, as a fraction of P and Kp.
%   'dP' (1e-6), 'dKp' (1e-5) minimum absolute perturbation sizes for P and
%      Kp to determine gradients. The actual perturbation is max(dPrel*P, dP)
%      and max(dKprel*Kp, dKp). Also 'dPabs' and 'dKpabs' are accepted.
%   'info' (1) how much detail to provide, 0 - 3
%   'criterion' (1e-6) criterion for stopping the gradient search
%      when internally looking for fixed points 
%   'parmvary1' ('P') and 'parmvary2' ('Kp') specify which variables to vary in
%      order to try to meet the target gait.  Any two parameters of a model
%      can be used, from get(w, 'parms'), but they should have separate
%      effects on speed and step length.  There are also options for
%      'dparm1rel', 'dparm2rel' (0.001), 'dparm1abs' and 'dparm2abs' (1e-5).
%   'relativedelta' (0) whether to use relative perturbations instead of abs
%   'NMax' (250) Max # of mesh points to use in bvp
%   'solver' (@bvp5c) Specify the handle of the bvp solver to use. Pre-
%          Matlab 2008b, use bvp4c instead of bvp5c
%
%   Additional output [w, cnvrg, sol] includes the solution structure, 
%   containing fields x, y, and parameters associated with the bvp.

% Changes
% Added by Art 6/2009 as a way to improve upon findgaitspeed

speed0 = []; steplength0 = []; stepfreq0 = []; opt_argin = {};
info = 1; stepsize = 1; parmvary = []; criterion = 1e-8; 
fgscriterion = criterion;
reldelta = 0; % Set to 1 if you want findgaitspeedb to always use relative delta perturbations
extrapolate = 1; 
reltol = 1e-4; NMax = 250; bvpsolver = @bvp5c;
maxrecursion = 8; recursion = 0;

if nargin == 0 || nargin == 2
 error('bvpgaitspeed(w) or bvpgaitspeed(w, speed, steplength)')
elseif nargin == 1
 speed0 = 0.4; % look for a nominal gait of equivalent to 1.25 m/s
 stepfreq0 = 1.8 / sqrt(9.81); % and 1.8 Hz
elseif nargin >= 3 && isnumeric([varargin{1} varargin{2}]) % numeric inputs for speed and step length
 speed0 = varargin{1}; steplength0 = varargin{2};
 if isempty(speed0), speed0 = 0.4; end % if either is empty, use default values
 if isempty(steplength0), steplength0 = speed0 /(1.8 / sqrt(9.81)); end
 opt_argin = varargin(3:end); % skip two arguments
else
 opt_argin = varargin;
end

walk = w0;

parmvary1 = 'P'; parmvary2 = 'Kp'; % usually vary P and Kp to get desired gait
dparm1rel = 0.001; dparm2rel = 0.001; % relative recommended changes
dparm1abs = 1e-6; dparm2abs = 1e-5; % minimum absolute changes

stepfreq0 = speed0 / steplength0; 

cnvrg = 0; caught = 0;

saveopt = opt_argin; % save for recursion
nopt = length(opt_argin); iopt = 1;
while iopt <= nopt-1, % make sure there is at least one pair left 
opt = opt_argin{iopt};
val = opt_argin{iopt+1};
  switch opt
    case 'speed'                           % findgaitspeed options
      speed0 = val;
    case 'steplength'
      steplength0 = val;
    case {'stepfrequency','stepfreq'}
      stepfreq0 = val;
    case 'parmvary1'
      parmvary1 = val;
    case 'parmvary2'
      parmvary2 = val;
    case {'dPrel','dparm1rel'}
      dparm1rel = val;
    case {'dKprel', 'dparm2rel'}
      dparm2rel = val;
    case {'dP','dparm1abs','dPabs'}  % minimum stepsize of findng gradient with dP
      dparm1abs = val;
    case {'dKp','dparm2abs','dKpabs'} % minimum stepsize of finding gradient with Kp
      dparm2abs = val;
    case 'relativedelta'
      reldelta = val;
    case 'stepsize' % stepsize of varying P, Kp
      stepsize = val; % criterion for ending findgaitspeed
    case 'fgscriterion'
      fgscriterion = val;
    case 'info'
      info = val;
    case 'delta'                       % findgait options
      delta = val;
    case 'criterion'                       % findgait options
      criterion = val;
    case 'extrapolate'
      extrapolate = val;
    case 'solver'
      bvpsolver = val;
    case 'NMax'
      NMax = val;
     case 'maxrecursion'
       maxrecursion = val;
    otherwise
      warning('bvp:options',['findgaitspeedb options: speed, steplength, dP, dKp, fgscriterion, stepfrequency, info,\n', ...
        'dPrel, dKprel, parmvary1, parmvary2, dparm1abs, dparm2abs, dparm1rel, dparm2rel\n',...
        'criterion, NMax, maxrecursion'])
  end % switch
   iopt = iopt + 2; % going through a pair at a time
end % while going through options

% Expect two of {speed0, steplength0, stepfreq0} to be specified, but
% fill in the third if one of them is empty
if isempty(speed0)
  speed0 = steplength0 * stepfreq0;
elseif isempty(steplength0)
  steplength0 = speed0 / stepfreq0;
elseif isempty(stepfreq0)
  stepfreq0 = speed0 / steplength0;
end

N = get(walk, 'N');
wstar = walk;

statoption = 'off';
if info >= 1
  fprintf(1,'findgaitspeedb target: speed = %g, step length = %g\n', speed0, steplength0);
  fprintf(1,'  parmvary1 = %s, parmvary2 = %s\n', parmvary1, parmvary2);
  if info > 2
    statoption = 'on';
  end
end

walkparms = get(walk, 'parms');

% we'll use the existing xstar as the fixed point.
x0 = get(walk, 'xstar');
if isnumeric(x0)
  solinit = bvpsetup(walk, x0); % must be provided for the walk class
elseif isstruct(x0) % looks like it is a structure suitable for solinit
  solinit = x0;
  x0 = solinit.y(:,1);
end

% augment the bvp parameters with the two parameters we've agreed to vary
solinit.parameters = [solinit.parameters walkparms.(parmvary1) walkparms.(parmvary2)];

options = bvpset('AbsTol', criterion, 'RelTol', reltol, 'Stats', statoption);

% check out the initial walk object to find its speed
[bcres, steplength, stepperiod] = bvpbcfun(walk, solinit);
currentgait = [steplength stepperiod];  % express targets as step length and period
targetgait = [steplength0 1/stepfreq0];

% This is where everything happens, in a recursive routine
[sol, cnvrg] = bvpspeedsearcher(currentgait, targetgait, solinit, recursion);

% evaluate whether we met the residual
wstar = set(walk, 'xstar', sol.y(:,1)', parmvary1, sol.parameters(end-1), parmvary2, sol.parameters(end));
[bcres, steplength, stepperiod] = bvpbcfun(wstar, sol);
cnvrg = all(abs(bcres) < criterion) && cnvrg;

if info > 0
  fprintf(1,'bvpsearchspeed return  cnvrg = %d  bcres = %g speed = %g, steplength = %g\n',...
    cnvrg, max(abs(bcres)), steplength/stepperiod, steplength);
  % bcres is the boundary condition residual, a measure of how well we
  % match all fixed point conditions x0(k+1) = x0, the Poincare section,
  % and the speed and step length.
end

%% bvpspeedsearcher (recursive routine)
function [sol, cnvrg] = bvpspeedsearcher(currentgait, targetgait, solinit, recursion)
% need to set up a search where we try to converge, and if it doesn't work
% we go a recursion level deeper
caught = 0;
steplength0 = targetgait(1); stepperiod0 = targetgait(2);
infoindent = repmat('  ', 1, recursion);

try
  sol = bvpsolver(@bvpspeedfwalk, @bvpspeedbcfun, solinit, options);
catch ME
%  if isempty(findstr('SingJac', ME.identifier)), % a real error (old code using deprecated findstr) 
  if isempty(strfind(ME.identifier, 'SingJac')), % a real error
    rethrow(ME); % anything but the singular Jacobian is just an error
  else % singular Jacobian, meaning poor convergence
    if info >=2
      fprintf('Singular Jacobian: ');
    end
    sol = solinit; % whereas we're relatively forgiving of a singular
    caught = 1;    % Jacobian, which is just like a failed search
  end
end % try

% evaluate bc's

newwalk = set(walk, parmvary1, sol.parameters(end-1), parmvary2, sol.parameters(end));
[bcres, steplength, stepperiod] = bvpbcfun(newwalk, sol);
cnvrg = all(abs(bcres) < criterion) && ~caught;
if info >=1 && ~cnvrg, fprintf('Initial bvp unsuccessful\n'); end;

if ~cnvrg,
  % so try bisection 1
  if info >= 2, fprintf('%sbisection 1: ', infoindent); end
  if recursion >= maxrecursion
    if info >= 2,
      fprintf('max recursion limit reached in bvpspeedsearcher, recursion = %d\n', recursion)
    end
    cnvrg = 0;
    sol = solinit;
    return;
  end
  % good to go a level deeper
  halftarget = mean([currentgait; targetgait]);
  [halfsol, halfcnvrg] = bvpspeedsearcher(currentgait, halftarget, solinit, recursion+1);
  if halfcnvrg == 1, % successful half-step, so try the next half: bisection 2
    if info >=2, fprintf('%sbisection 2: ', infoindent); end;
    [sol, cnvrg] = bvpspeedsearcher(halftarget, targetgait, halfsol, recursion+1);
  else % unsuccessful half-step, so quit out with cnvrg = 0
    cnvrg = halfcnvrg;
    sol = solinit;
  end
end % if ~cnvrg
if cnvrg
  fprintf(1,'bvp successful search, recursion level = %d\n', recursion);
else % unsuccessful, quit out
  fprintf(1,'bvp unsuccessful, recursion level = %d\n', recursion);
end %
      
  %% bvpspeedfwalk nested
  function xdot = bvpspeedfwalk(varargin)
    % BVPSPEEDFWALK automatically calls fwalk in such a way that the time is normalized
    % and the parameters are passed appropriately.
    % Typical call: bvpspeedfwalk(t, x, unknownbvpparms) or
    % bvpspeedfwalk(t, x, region, unknownbvpparms) for multipoint bvps
    
    % the last entry in varargin is assumed to be the unknown parms:
    unknownbvpparms = varargin{end};
    
    % we'll send fwalk our own set of slightly different parameters
    passedparms = walkparms;
    passedparms.(parmvary1) = unknownbvpparms(end-1);
    passedparms.(parmvary2) = unknownbvpparms(end);
    
    % and we'll just call bvpfwalk with the same arguments, except
    % augmenting them with the walk object and the passed parameters
    xdot = bvpfwalk(varargin{:}, walk, passedparms);
  end % bvpspeedfwalk

  %% bvpspeedbcfun nested
  function bcres = bvpspeedbcfun(varargin)
  % BVPSPEEDBCFUN automatically calls bvpbcfun in such a way that the time
  % is normalized and the appropriate parameters are passed. It calls the
  % regular bvpbcfun to find the fixed point boundary conditions, and 
  % augments them with extra conditions for desired step length and period.
  % It also extracts the unknown parameters of the bvp and applies them
  % to the "known" parameters of the walk object, so that they are used
  % in the equations of motion.
  % Typical call:  bvpspeedbcfun(xa, xb, unknownparms), with nested
  % variables steplength0, stepperiod0.
  % Alternative call: bvpspeedbcfun([any arguments], unknownparms)
  % where [any arguments] is passed straight to bvpbcfun, which accepts
  % more than one method of calling.
    
    % the last entry in varargin is assumed to be the unknown bvp parms:
    unknownbvpparms = varargin{end};
    % these contain the regular bvpsearch parameters (such as step period
    % T) and are augmented with two others for parmvary1 and parmvary2
    % which will be the last two entries of the bvp parameters array.
    
    % we'll send bcfun our own set of slightly different walk parameters
    passedparms = walkparms;
    passedparms.(parmvary1) = unknownbvpparms(end-1);
    passedparms.(parmvary2) = unknownbvpparms(end);
    
    % and we'll just call bvpfwalk with the same arguments, except
    % augmenting them with the walk object and the passed parameters
    [bcres, steplength, stepperiod] = bvpbcfun(varargin{:}, walk, passedparms);
    % we also augment the boundary conditions with the step length and
    % period conditions
    bcres = [bcres; steplength - steplength0; stepperiod - stepperiod0];
    
  end % bvpspeedbcfun
      
end % bvpspeedsearcher


end % findgaitspeedb