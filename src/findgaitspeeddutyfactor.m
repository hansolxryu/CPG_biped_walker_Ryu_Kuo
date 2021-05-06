function [w, cnvrg] = findgaitspeeddutyfactor(w0, varargin)
% FINDGAITSPEEDDUTYFACTOR determines the model parameters needed
%   to produce a desired speed and steplength and dutyfactor.
%   [w, cnvrg] = findgaitspeeddutyfactor(w0, 'speed', speed, 'steplength', steplength, 'dutyfactor', dutyfactor)
%   takes the existing gait w0, and tries to find parameter values (E,Kp,Kleg by default) to match
%   the desired speed, steplength and duty factor, returning the new gait in w, with the
%   flag cnvrg indicating whether the search was successful
%   [w, cnvrg] = findgaitspeeddutyfactor(w0, speed, steplength, dutyfactor) takes the known
%   existing gait and tries to find parameter values (E,Kp,Kleg by default) to match the desired
%   speed, steplength and duty factor, returning the new gait in w, with the flag
%   cnvrg indicating whether the search was successful.
%   findgaitspeeddutyfactor(w0, speed, steplength, dutyfact, option, val, ...) applies the
%   desired options, which include the following:
%   'dparm1rel', 'dparm2rel', 'dparm3rel' (0.001) relative perturbation sizes parameters to
%      determine gradients.
%   'dparm1abs' (1e-6), 'dparm2abs' (1e-5), 'dparm3abs' (1e-5) minimum absolute perturbation sizes 
%      parameters to determine gradients. The actual perturbation is
%      max(dparmrel*parmvalue, dparmabs)
%   'info' (2) how much detail to provide, 0 - 3
%   'fgscriterion' (1e-4) how closely speed, steplength and duty factor should match
%   'criterion' (1e-7) criterion for stopping the gradient search
%      when internally looking for fixed points (gradsearch option)
%   'stepsize' (0.5) how quickly to step the adjustments in parameters
%      relative to the Newton's method change
%   'extrapolate' (1) whether to try to extrapolate new fixed points
%      based on previous data, generally supposed to be faster but not
%      completely reliable
%   'parmvary1' ('E') and 'parmvary2' ('Kp') and 'parmvary3' ('Kleg') specify which variables to vary in
%      order to try to meet the target gait.  Any three parameters of a model
%      can be used, from get(w, 'parms'), but they should have separate
%      effects on speed and step length.  

% Shawn O'Connor

speed0 = []; steplength0 = []; stepfreq0 = []; dutyfact0 = []; opt_argin = {};

if nargin == 0 || nargin == 2
 error('findgaitspeed(w) or findgaitspeed(w, speed, steplength, dstime)')
elseif nargin == 1
 speed0 = 0.4; % look for a nominal gait of equivalent to 1.25 m/s
 stepfreq0 = 1.8 / sqrt(9.81); % and 1.8 Hz
 dutyfact0 = 0.6;
elseif nargin >= 4 && isnumeric([varargin{1} varargin{2} varargin{3}]) % numeric inputs for speed and step length
 speed0 = varargin{1}; steplength0 = varargin{2}; dutyfact0 = varargin{3};
 if isempty(speed0), speed0 = 0.4; end % if either is empty, use default values
 if isempty(steplength0), steplength0 = speed0 /(1.8 / sqrt(9.81)); end
 if isempty(dutyfact0), dutyfact0 = 0.6; end % if either is empty, use default values
 opt_argin = varargin(4:end); % skip two arguments
else
 opt_argin = varargin;
end

w = w0;
stepsize = 0.50;
reldelta = 0; % Set to 1 if you want findgaitspeed to always use relative delta perturbations
criterion = 1e-10; extrapolate = 1; fgscriterion = 1e-3;
info = 2;
parmvary1 = 'E';    parmvary2 = 'Kp';   parmvary3 = 'Kleg'; % usually vary P and Kp to get desired gait
dparm1rel = 0.0005; dparm2rel = 0.001;  dparm3rel = 0.001; % relative recommended changes
dparm1abs = 1e-6;   dparm2abs = 1e-5;   dparm3abs = 1e-5; % minimum absolute changes

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
   case 'dutyfact'
     dutyfact0 = val;
   case 'parmsearch'
     findgaitoptions = {findgaitoptions{:}, opt, val};
   case 'parmvary1'
     parmvary1 = val;
   case 'parmvary2'
     parmvary2 = val;
   case 'parmvary3'
     parmvary3 = val;
   case 'dparm1rel'
     dparm1rel = val;
   case 'dparm2rel'
     dparm2rel = val;
   case 'dparm3rel'
     dparm3rel = val;
   case 'dparm1abs'
     dparm1abs = val;
   case 'dparm2abs'
     dparm2abs = val;
   case 'dparm3abs'
     dparm3abs = val;
   case 'relativedelta'
     reldelta = val;     
   case 'stepsize' % stepsize of varying parameters
     stepsize = val; % criterion for ending findgaitspeed  
   case 'fgscriterion'
     fgscriterion = val;
   case 'info'
     info = val;
   case 'criterion'                       % findgait options
     criterion = val;
   case 'extrapolate'
     extrapolate = val;
   case 'parmgrad'     
   otherwise
     warning(['findgaitspeed options: speed, steplength, dutyfact, stepsize, fgscriterion, stepfrequency, info,\n', ...
         'parmvary1, parmvary2, parmvary3, dparm1abs, dparm2abs, dparm3abs, dparm1rel, dparm2rel, dparm3rel\n',...
         'findgait options: criterion, extrapolate'])
   end
end

if isempty([speed0 steplength0 stepfreq0 dutyfact0])
 speed0 = 0.40; stepfreq0 = 1.8 /sqrt(9.81);
 steplength0 = speed0 /stepfreq0;
 dutyfact0 = 0.6;
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
elseif isempty(dutyfact0)
 error('findgaitspeed requires duty factor to be set')
end

ginfo = info - 1;
enew = 0; e = 1;  cnvrg = 0;
i = 1;
[v, sl, sf, df] = gaitspeed(w);
while max(abs(e)) > fgscriterion
 parm1 = get(w, parmvary1); parm2 = get(w, parmvary2); parm3 = get(w, parmvary3);
 if reldelta
     dparm1 = dparm1rel*parm1; dparm2 = dparm2rel*parm2; dparm3 = dparm3rel*parm3;
 else
     dparm1 = max(dparm1rel*parm1, dparm1abs); dparm2 = max(dparm2rel*parm2, dparm2abs); dparm3 = max(dparm3rel*parm3, dparm3abs);
 end
 % Find gradient of step length, step frequency, and duty factor as a function of the three parameters
 w1 = set(w, parmvary1, parm1+dparm1); w2 = set(w, parmvary2, parm2+dparm2);  w3 = set(w, parmvary3, parm3+dparm3);
 [w1, cnvrg1] = findgait(w, w1, [], 'criterion', criterion, 'parmvary', parmvary1, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});
 [w2, cnvrg2] = findgait(w, w2, [], 'criterion', criterion, 'parmvary', parmvary2, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});
 [w3, cnvrg3] = findgait(w, w3, [], 'criterion', criterion, 'parmvary', parmvary3, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});
 if ~cnvrg1
   fprintf(1, '  could not converge on delta %s\n', parmvary1); cnvrg = 0; break;
 elseif ~cnvrg2
   fprintf(1, '  could not converge on delta %s\n', parmvary2); cnvrg = 0; break;
 elseif ~cnvrg3
   fprintf(1, '  could not converge on delta %s\n', parmvary3); cnvrg = 0; break;   
 end
 [v1, sl1, sf1, df1] = gaitspeed(w1); [v2, sl2, sf2, df2] = gaitspeed(w2);  [v3, sl3, sf3, df3] = gaitspeed(w3);
 dy = [(sl1-sl) (sl2-sl) (sl3-sl); (sf1-sf) (sf2-sf) (sf3-sf); (df1-df) (df2-df) (df3-df)];
 J = dy ./ [dparm1 dparm2 dparm3; dparm1 dparm2 dparm3; dparm1 dparm2 dparm3]; % gradient of gait parameters wrt model parms
 % Adjust variable parameters to make a new guess
 deltaparms = -stepsize * (J \ [sl-steplength0; sf-stepfreq0; df-dutyfact0]);
 wnew = set(w, parmvary1, parm1+deltaparms(1), parmvary2, parm2+deltaparms(2), parmvary3, parm3+deltaparms(3));
 maxdeltaparms = max(abs(deltaparms));   
 if maxdeltaparms == abs(deltaparms(1))
   [wnew,cnvrg] = findgait(w, wnew, [], 'parmvary', parmvary1, 'criterion', criterion, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});
 elseif maxdeltaparms == abs(deltaparms(2))
   [wnew,cnvrg] = findgait(w, wnew, [], 'parmvary', parmvary2, 'criterion', criterion, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});
 elseif maxdeltaparms == abs(deltaparms(3))
   [wnew,cnvrg] = findgait(w, wnew, [], 'parmvary', parmvary3, 'criterion', criterion, 'extrapolate', extrapolate, 'info', ginfo, findgaitoptions{:});  
 end
 if ~cnvrg
   fprintf(1,'  could not converge on new guess\n');
   break
 end
 [v, sl, sf, df] = gaitspeed(wnew);
 enew = [sl-steplength0 sf-stepfreq0 df-dutyfact0];
 
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
fprintf(1, 'findgaitspeed took %d steps ', i);
if cnvrg
 fprintf(1, 'success\n');
end