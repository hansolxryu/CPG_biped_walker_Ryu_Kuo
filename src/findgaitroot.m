function [wstar,cnvrg,errmsg,wsuccess] = findgaitroot(w1, w2, varargin)
% FINDGAITROOT Finds a fixed point for a walking model by doing
% a creeping gradient search, starting from a known gait.
% This version differs from findgait in that it calls rootsearch,
% which allows a user-specified rootfunction.
%
% [wstar, cnvrg, errmsg, wsuccess] = findgaitroot(w1, w2, x0) takes the model w2, and returns the
% fixed point in model wstar, which includes xstar.  The search starts
% with the known gait w1, and optionally an initial guess x0.
% cnvrg indicates whether the search worked or not
% errmsg provides the last error encountered
% if the search fails, it returns the last successful gait closest to w21
% [wstar, cnvrg,errmsg, wsuccess] = findgait... returns an array of gaits
% that includes all the successful intermediate searches made between
% w1 and w2.
% 'parmvary' ('gamma') is the parameter being varied, which helps
%    both for printing progress info, and especially if extrapolation
%    is used
% 'info' (1) determines whether progress is printed to the screen
% 'extrapolate' (1) tells whether to save intermediate gaits and use them
%    to extrapolate and find the new gaits using linear or quadratic fits
% 'minparmdiff'
% Other options may be included for rootsearch, which is called by findgait.
% 
% [...] = findgait(w1, w2) or findgait(w1, w2, 'parmvary', ...) also works
%   In other words, the initial guess can be omitted entirely.
%
% rootsearch options:
% 'delta' (1e-5) is the perturbation size for evaluating gradients
% 'criterion' (1e-10) is the criterion for stopping the search
% 'stepsize' (1) is the size of steps in Newton's method, which can
%            be scaled to less than one to slow the search
%  
% 'rootfunction' the function whose root is to be found
% 'roottarget'   the value that constitutes a target for the root

% This function requires interpolategait and extrapolategait

% Art Kuo

% Changes
%   v 1.1 added minparmdiff parameter to avoid unending
%         interpolations (05/07/2004)
%   v 1.2 added error message reporting last error
%         wstar now returns the last successful gait if there is no
%         convergence (07/18/2007)
%   v 1.3 added default delta and criterion detection 
%         added non-numeric parameter handling (07/04/2008)
%         
%   v 1.4 has several changes. First, in the event of an unsuccessful
%         search, wstar now returns the final successful gait (the last
%         one in wsuccess, not the first or w1). Second, findgait can 
%         now be called without an initial guess argument, and instead
%         that entire argument can be skipped (no need to use []). 
%         Finally, there is a hack for dealing with hybridsys classes
%         that has been forked. The top-level findgait no longer cares
%         about hybrid systems, and a more specialized version is 
%         included in hybridsys. (modified by Art 07/23/2008)

% Example: findgait(walk2('normal'),set(walk2('normal'),'gamma',1))

% Things to do: errmsg is not particularly useful, as it only captures
% errors from rootsearch. It would help to return more informative
% information

drawnow;
cnvrg = 0; wsuccess = []; errmsg = ''; extrapolate = 1; 

delta = 1e-5; info = 1; criterion = 1e-10; stepsize = 1; parmvary = []; minparmdiff = 1e-5;

% third argument (first varargin) could be x0, or else it could be first
% of a list of option pairs.
if nargin >= 3 && isnumeric(varargin{1})
  x0 = varargin{1};
  varargin = varargin(2:end);
elseif nargin < 3 || isempty(varargin{1}) || ~isnumeric(varargin{1})  % x0 not given, get it from xstar
 x0 = get(w1, 'xstar');
 if isempty(x0)
   error('findgait requires an initial guess x0')
 end
end

% Deal with the parameter list
options = varargin;
opt_argin = options;
rootsearchoptions = {};
while length(opt_argin) >= 2,
 opt = opt_argin{1};
 val = opt_argin{2};
 opt_argin = opt_argin(3:end);
 switch opt
   case {'delta','criterion','stepsize'}
     % options for rootsearch, do nothing
   case 'info'
     info = val;
   case 'parmvary'
     parmvary = val; % also need to tell rootsearch which parm to vary
     rootsearchoptions = {rootsearchoptions{:}, opt, val};
   case 'parmignore'
     parmignore = val;
   case 'minparmdiff'
     minparmdiff = val;
   case 'extrapolate'
     extrapolate = 1;
     if ischar(extrapolate) % user entered a string
       parmvary = val;      % which is the parameter being varied
     end
   case {'AbsTol','RelTol', 'tf','rootfunction', 'roottarget', 'parmvary1', 'parmvary2'}
     rootsearchoptions = {rootsearchoptions{:}, opt, val};
     % options for rootsearch, do nothing
   otherwise
     warning('Findgaitroot options: extrapolate, parmvary, delta, info, criterion, stepsize')
 end
end

% If the user didn't specify which parameter is being varied, pick the
% parameter which undergoes the largest change
if isempty(parmvary)
 parms2 = get(w2, 'parms'); fnames = fieldnames(parms2);
 % allow for non-numeric parameters (which cannot be changed by findgait)
 pv1 = struct2cell(get(w1,'parms'));  pv2 = struct2cell(get(w2,'parms'));  i = 1; 
 while i <= length(pv1); 
     if ~isnumeric(pv1{i});  pv1 = cat(1,pv1(1:i-1),pv1(i+1:end));  pv2 = cat(1,pv2(1:i-1),pv2(i+1:end));  fnames = cat(1,fnames(1:i-1),fnames(i+1:end));  end;  
     i = i+1;
 end
 pvalues1 = cell2mat(pv1); % a numeric array of all parameter values
 pvalues2 = cell2mat(pv2);
 pdiff = abs(pvalues2 - pvalues1);
 plarge = max(abs(pvalues2),abs(pvalues1)); plarge(plarge==0) = Inf; % ignore 0 parameters
 [p, iparm] = max(pdiff ./ plarge);
 if isempty(iparm)
   error('Findgait could not determine which parameter is being varied')
 end
 parmvary = fnames{iparm};
 options = {options{:}, 'parmvary', parmvary}; % add this knowledge to parameters
else % we know exactly which parm is being varied
 pvalues1 = get(w1, parmvary); pvalues2 = get(w2, parmvary);
%  p = abs(pvalues1-pvalues2)/max(abs(pvalues1),abs(pvalues2));
 p = abs(pvalues1-pvalues2)/max([abs(pvalues1),abs(pvalues2),any(pvalues1==0),any(pvalues2==0)]);
end

% start with the gait w1 and see if its fixed point works for w2

if info > 0  % diagnostic information
 fprintf(1,'findgait start: ');
 if isempty(parmvary)
   display(w1);
 else
   fprintf(1, '%s = %g\n', parmvary, get(w1,parmvary));
 end
 fprintf(1,'findgait target: ');
 if isempty(parmvary)
   display(w2);
 else
   fprintf(1, '%s = %g\n', parmvary, get(w2,parmvary));
 end
end

% The first thing to try is to converge on the target directly
try
 [wstar, cnvrg] = rootsearch(w2, x0, rootsearchoptions{:});
catch
 errmsg = lasterr;
 fprintf(lasterr);
 if info > 0
   fprintf(1,'rootsearch 1 failed\n')
 end
 % if we had trouble converging AND there's basically no difference
 % between the current and the target gait, time to get out
 if p < minparmdiff % there's basically no difference with the current gait
   % and we shouldn't bother with interpolating further
   fprintf(1, 'interpolation limit reached: minparmdiff\n');
   wstar = w1; cnvrg = 0;
   return
 end
end

if cnvrg, % if you're already done, return successfully
 if info > 0
   fprintf('rootsearch 1 succeed, findgait return\n');
 end
 %  wsuccess = wstar;
 drawnow;
 return
else
   if info > 0
     fprintf(1, 'parmdiff = %g\n',p);
   end
 if p < minparmdiff % there's basically no difference with the current gait
   % and we shouldn't bother with interpolating further
   fprintf(1, 'interpolation limit reached: minparmdiff\n');
   wstar = w1; cnvrg = 0;
   return
 end   
end

% If direct convergence doesn't work, then it's time to interpolate
% w1 works, w2 doesn't, so let's interpolate
if info > 0, fprintf(1,'interpolate begin\n'); end
[wintermediate, cnvrg, errmsg, wsuccess] = findgaitroot(w1, interpolategait(w1,w2),[],options{:});
if info > 0, fprintf(1,'interpolate end\n'); end

if cnvrg, % the interpolation succeeded, so try to get the target again
%           this time starting from the intermediate solution
  wsuccess = [wintermediate; wsuccess];
  cnvrg = 0;
  if extrapolate && length(wsuccess) >= 1 % try extrapolating
    xnew = extrapolategait([wsuccess; w1], w2, 'parmvary', parmvary);
    if info > 0, fprintf(1,'  extrapolating x = '); disp(xnew);  end;
  else
    xnew = [];
  end
  [wstar, cnvrg, errmsg, wothers] = findgaitroot(wintermediate, w2, xnew, options{:});
  wsuccess = [wsuccess; wothers];
else
 if info > 0, fprintf(1,'no convergence on intermediate target\n'); end
 %wstar = w1;
 wstar = wintermediate; % try to return the one that converged
end

if info > 0, fprintf('exit findgaitroot\n'); end