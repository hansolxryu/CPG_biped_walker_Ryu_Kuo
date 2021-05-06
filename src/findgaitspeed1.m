function [w, cnvrg] = findgaitspeed1(w0, varargin)
% FINDGAITSPEED1 determines the impulsive-pushoff or hip spring needed
%   to produce a desired speed or steplength.
%   [w, cnvrg] = findgaitspeed1(w0, 'speed', speed)
%   takes the existing gait w0, and tries to find P values to match
%   the desired speed, returning the new gait in w, with the
%   flag cnvrg indicating whether the search was successful.
%   [w, cnvrg] = findgaitspeed(w0, 'steplength', sl) or 'stepfreq', sf
%   uses different targets for the gait.
%   [w, cnvrg] = findgaitspeed1(w0) automatically assumes speed is the
%   target, and uses nominal v = 0.4.
%   [w, cnvrg] = findgaitspeed1(w0, speed) takes the known
%   existing gait and tries to find P and Kp values to match the desired
%   speed, returning the new gait in w, with the flag
%   cnvrg indicating whether the search was successful.
%   findgaitspeed1(w0, speed, option, val, ...) applies the
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
%   'parmvary' ('P') specify which variable (such as Kp) to vary in
%      order to try to meet the target gait.  Any single parameter of a model
%      can be used, from get(w, 'parms'), but it should have an effect on
%      the desired target value of speed. There are also options for
%      'dparmrel',  (0.001), 'dparmabs' 1e-5) for relative or absolute
%      perturbations to parameter.
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
%   v 1.3b copied findgaitspeed into findgaitspeed1, for just one
%     parameter (Art 8/2011)

opt_argin = {};
target = 'speed'; % could also be 'steplength' or 'stepfreq' as returned
% as fields in the output of walkspeed(w).
targetval = 0.4; % default speed to aim for, equiv to 1.25 m/s human

if nargin == 0
 error('findgaitspeed(w) or findgaitspeed(w, speed)')
elseif nargin == 1 % just given a gait, find the nominal gait
  target = 'speed';
  targetval = 0.4; % look for a nominal gait of equivalent to 1.25 m/s
 %stepfreq0 = 1.8 / sqrt(9.81); % and 1.8 Hz
elseif nargin >= 2 && isnumeric(varargin{1}) % numeric input for speed
 targetval = varargin{1}; 
 opt_argin = varargin(2:end); % skip an argument
else
 opt_argin = varargin;
end

w = w0;
stepsize = 0.5; 
reldelta = 0; % Set to 1 if you want findgaitspeed to always use relative delta perturbations
delta = 1e-8; criterion = 1e-9; extrapolate = 1; fgscriterion = 1e-3;
info = 2;
parmvary = 'P'; parmvary2 = 'Kp'; % usually vary P and Kp to get desired gait
dparmrel = 0.001;  % relative recommended changes
dparmabs = 1e-6;  % minimum absolute changes

while length(opt_argin) >= 2,
 opt = opt_argin{1};
 val = opt_argin{2};
 opt_argin = opt_argin(3:end);
 switch opt
   case {'speed', 'steplength', 'stepfreq'}        % findgaitspeed options
     target = opt;
     targetval = val;
   case 'parmvary'
     parmvary = val;
   case {'dPrel','dparmrel'}
     dparmrel = val;
   case {'dP','dparmabs','dPabs'}  % minimum stepsize of findng gradient with dP
     dparmabs = val;
   case 'relativedelta'
     reldelta = val;
   case 'stepsize' % stepsize of varying P
     stepsize = val; % criterion for ending findgaitspeed1
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
     warning(['findgaitspeed1 options: speed, steplength, stepfreq, dP, dKp, stepsize, fgscriterion, stepfrequency, info,\n', ...
         'dPrel, parmvary, dparmabs, dparmrel, \n',...
         'findgait options: delta, criterion, extrapolate'])
   end
end


if isempty(target) || isempty(targetval)
 error('findgaitspeed requires one of speed, steplength, stepfreq to be set')
end

ginfo = info - 1;
enew = 0; e = 1;  cnvrg = 0;
i = 1;
speedinfo = walkspeed(w); % speedinfo will have fields 'speed', 'steplength', 'stepfreq'
while max(abs(e)) > fgscriterion
 parm = get(w, parmvary); 
 if reldelta
   dparm = dparmrel*parm; 
 else
   dparm = max(dparmrel*parm, dparmabs); 
 end
 % Find gradient of target as a function of
 % the parameter (P by default)
 w1 = set(w, parmvary, parm+dparm);
 [w1, cnvrg1] = findgait(w, w1, [], 'criterion', criterion, 'parmvary', parmvary, 'extrapolate', extrapolate, 'info', ginfo);
 if ~cnvrg1
   fprintf(1, '  could not converge on delta %s\n', parmvary); cnvrg = 0; break;
 end
 speedinfo1 = walkspeed(w1);
 dy = speedinfo1.(target) - speedinfo.(target);
 J = dy ./ dparm; % gradient of target value wrt parm
 % Adjust variable parameters to make a new guess
 deltaparm = -stepsize * (J \ (speedinfo.(target)-targetval)); 
 wnew = set(w, parmvary, parm+deltaparm);
 [wnew,cnvrg] = findgait(w, wnew, [], 'parmvary', parmvary, 'criterion', criterion, 'extrapolate', extrapolate, 'info', ginfo);
 if ~cnvrg
   fprintf(1,'  could not converge on new guess\n');
   break
 end
 speedinfo = walkspeed(wnew);
 enew = speedinfo.(target)-targetval;
 if max(abs(enew)) > max(abs(e))
   cnvrg = 0;
   if info > 0
     fprintf(1,' findgaitspeed1 could not reduce error on gait\n');
   end
   break
 end
 cnvrg = 1;
 e = enew;
 w = wnew; i = i + 1;
 if info
   fprintf(1, '  e = %g\n', max(abs(e)));
 end
end % while

if info
  fprintf(1, 'findgaitspeed1 took %d steps ', i);
  if cnvrg
    fprintf(1, 'success\n');
  end
end % if info
if ~cnvrg
  warning('findgaitspeed1 did not converge')
end

end % function findgaitspeed1



