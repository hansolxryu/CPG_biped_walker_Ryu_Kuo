function [w, cnvrg] = findgaitspeed(w0, varargin)
% FINDGAITSPEED determines the impulsive-pushoff and hip spring needed
%   to produce a desired speed and steplength.
%   [w, cnvrg] = findgaitspeed(w0, 'speed', speed, 'steplength', steplength)
%   takes the existing gait w0, and tries to find P and Kp values to match
%   the desired speed and steplength, returning the new gait in w, with the
%   flag cnvrg indicating whether the search was successful
%   [w, cnvrg] = findgaitspeed(w0, speed, steplength) takes the known
%   existing gait and tries to find P and Kp values to match the desired
%   speed and steplength, returning the new gait in w, with the flag
%   cnvrg indicating whether the search was successful.
%   findgaitspeed(w0, speed, steplength, option, val, ...) applies the
%   desired options, which include the following:
%   'dPrel', 'dKprel' (0.001) relative perturbation sizes for P and Kp to
%      determine gradients, as a fraction of P and Kp.
%   'dP' (1e-6), 'dKp' (1e-5) minimum absolute perturbation sizes for P and
%      Kp to determine gradients. The actual perturbation is max(dPrel*P, dP)
%      and max(dKprel*Kp, dKp). Also 'dPabs' and 'dKpabs' are accepted.
%   'info' (2) how much detail to provide, 0 - 3
%   'fgscriterion' (1e-4) how closely speed and steplength should match
%   'criterion' (1e-7) criterion for stopping the gradient search
%      when internally looking for fixed points (gradsearch option)
%   'stepsize' (0.5) how quickly to step the adjustments in P and Kp
%      relative to the Newton's method change
%   'extrapolate' (1) whether to try to extrapolate new fixed points
%      based on previous data, generally supposed to be faster but not
%      completely reliable
%   'parmvary1' ('P') and 'parmvary2' ('Kp') specify which variables to vary in
%      order to try to meet the target gait.  Any two parameters of a model
%      can be used, from get(w, 'parms'), but they should have separate
%      effects on speed and step length.  There are also options for
%      'dparm1rel', 'dparm2rel' (0.001), 'dparm1abs' and 'dparm2abs' (1e-5).
%   'relativedelta' (0) whether to use relative perturbations instead of abs

% Art Kuo

% Changes
%   v 1.1 fixed bug in step freq, added 'speed', 'step
%     frequency', 'step length' as parameters (Art Kuo, 04/07/2004)
%   v 1.1a fixed bug for calling without any speed or step frequency
%     (04/08/2004)
%   v 1.2 added facility for percentage changes, as well as varying
%     any two parameters instead of just P and Kp (05/07/2004)
%   v 1.3 added default delta, criterion, and parmvary1 detection, 
%     as well as delta as a possible varargin (07/04/2008)
%   v 1.3a included relativedelta in help comment; note delta is 
%     probably not used correctly. (Art 07/24/2008)

speed0 = []; steplength0 = []; stepfreq0 = []; opt_argin = {};

if nargin == 0 || nargin == 2
 error('findgaitspeed(w) or findgaitspeed(w, speed, steplength)')
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

w = w0;
stepsize = 0.5; 
reldelta = 0; % Set to 1 if you want findgaitspeed to always use relative delta perturbations
delta = 1e-8; criterion = 1e-9; extrapolate = 1; fgscriterion = 1e-3;
info = 2;
parmvary1 = 'P'; parmvary2 = 'Kp'; % usually vary P and Kp to get desired gait
dparm1rel = 0.001; dparm2rel = 0.001; % relative recommended changes
dparm1abs = 1e-6; dparm2abs = 1e-5; % minimum absolute changes
% If walk object has a property "ddelta" for the default delta size, e.g. because of the limitations of ADAMS, use it.  Same for criterion and parmvary1
try delta = get(w,'ddelta'); end;  try criterion = get(w,'dcrit'); end;  try parmvary1 = get(w,'dparmvary1'); end

findgaitoptions = {};
while length(opt_argin) >= 2,
 opt = opt_argin{1};
 val = opt_argin{2};
 opt_argin = opt_argin(3:end);
 switch opt
   case 'speed'                           % findgaitspeed options
     speed0 = val;
   case 'steplength'
     steplength0 = val;
   case {'stepfrequency','stepfreq'}
     stepfreq0 = val;
   case 'parmsearch'
     findgaitoptions = {findgaitoptions{:}, opt, val};
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
   otherwise
     warning(['findgaitspeed options: speed, steplength, dP, dKp, stepsize, fgscriterion, stepfrequency, info,\n', ...
         'dPrel, dKprel, parmvary1, parmvary2, dparm1abs, dparm2abs, dparm1rel, dparm2rel\n',...
         'findgait options: delta, criterion, extrapolate'])
   end
end

if isempty([speed0 steplength0 stepfreq0])
 speed0 = 0.40; stepfreq0 = 1.8 /sqrt(9.81);
 steplength0 = speed0 /stepfreq0;
end

if isempty(speed0)
 speed0 = steplength0 * stepfreq0;
elseif isempty(steplength0)
 steplength0 = speed0 / stepfreq0;
elseif isempty(stepfreq0)
 stepfreq0 = speed0 / steplength0;
end

stepfreq0 = speed0 / steplength0;

if isempty(speed0) || isempty(stepfreq0) || isempty(steplength0)
 error('findgaitspeed requires two of speed, steplength, stepfrequency to be set')
end

ginfo = info - 1;
enew = 0; e = 1;  cnvrg = 0;
i = 1;
[v, sl, sf] = gaitspeed(w);
while max(abs(e)) > fgscriterion
 parm1 = get(w, parmvary1); parm2 = get(w, parmvary2);
 if reldelta
     dparm1 = dparm1rel*parm1; dparm2 = dparm2rel*parm2;
 else
     dparm1 = max(dparm1rel*parm1, dparm1abs); dparm2 = max(dparm2rel*parm2, dparm2abs);
 end
 % Find gradient of step length and step frequency as a function of
 % the two parameters (P and Kp by default)
 w1 = set(w, parmvary1, parm1+dparm1); w2 = set(w, parmvary2, parm2+dparm2);
    [w1, cnvrg1] = findgait(w, w1, [], 'criterion', criterion, 'parmvary', parmvary1, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});
    [w2, cnvrg2] = findgait(w, w2, [], 'criterion', criterion, 'parmvary', parmvary2, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});

 if ~cnvrg1
   fprintf(1, '  could not converge on delta %s\n', parmvary1); cnvrg = 0; break;
 end
 if ~cnvrg2
   fprintf(1, '  could not converge on delta %s\n', parmvary2); cnvrg = 0; break;
 end
 [v1, sl1, sf1] = gaitspeed(w1);     [v2, sl2, sf2] = gaitspeed(w2);
 dy = [sl1-sl (sl2-sl); sf1-sf (sf2-sf)];
 J = dy ./ [dparm1 dparm2; dparm1 dparm2]; % gradient of step length and freq wrt parms
 % Adjust variable parameters to make a new guess
 deltaparms = -stepsize * (J \ [sl-steplength0; sf-stepfreq0]);
 wnew = set(w, parmvary1, parm1+deltaparms(1), parmvary2, parm2+deltaparms(2));
 if abs(deltaparms(1)) > abs(deltaparms(2))
   [wnew,cnvrg] = findgait(w, wnew, [], 'parmvary', parmvary1, 'criterion', criterion, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});
 else
   [wnew,cnvrg] = findgait(w, wnew, [], 'parmvary', parmvary2, 'criterion', criterion, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});
 end
 if ~cnvrg
   fprintf(1,'  could not converge on new guess\n');
   break
 end
 [v, sl, sf] = gaitspeed(wnew);
 enew = [sl-steplength0 sf-stepfreq0];
 if max(abs(enew)) > max(abs(e))
   cnvrg = 0;
   if info > 0
     fprintf(1,' findgaitspeed could not reduce error on gait\n');
   end
   break
 end
 cnvrg = 1;
 e = enew;
 w = wnew; i = i + 1;
 if info
   fprintf(1, '  e = %g\n', max(abs(e)));
 end
end
if info > 0
  fprintf(1, 'findgaitspeed took %d steps ', i);
  if cnvrg
    fprintf(1, 'success\n');
  end
end
