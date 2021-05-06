function [wstar,cnvrg,sol] = findgaitb(w1, w2, varargin)
% FINDGAITB Finds a fixed point for a walking model by doing
% a creeping bvp search, starting from a known gait.
% [wstar, cnvrg, sol] = findgaitb(w1, w2, x0) takes the model w2, and returns the
% fixed point in model wstar, which includes xstar.  The search starts
% with the known gait w1, and optionally an initial guess x0.
% cnvrg indicates whether the search worked or not
% errmsg provides the last error encountered
% if the search fails, it returns the last successful gait closest to w2
% 'parmvary' ('gamma') is the parameter being varied, which helps
%    both for printing progress info, and especially if extrapolation
%    is used
% 'info' (1) determines whether progress is printed to the screen
% 'maxrecursion' (8) determines maximum number of levels of recursion
% 'minparmdiff'
% Other options may be included for bvpsearch, which is called by findgait.
% 
% [...] = findgait(w1, w2) or findgait(w1, w2, 'parmvary', ...) also works
%   In other words, the initial guess can be omitted entirely.
%
% bvpsearch options:
% 'criterion' (1e-8) is the criterion for stopping the search
% 'NMax'      (250) Max # of mesh points to use in bvp
% 'solver' (@bvp5c2) Specify the handle of the bvp solver to use. Pre-
%          Matlab 2008b, use bvp4c instead of bvp5c

% This function requires interpolategait and extrapolategait

% Art Kuo

% Changes
%   Art added findgaitb(vp), based on findgait but for bvps

% Example: findgaitb(walk2('normal'),set(walk2('normal'),'gamma',1))

drawnow;
cnvrg = 0;  
recursion = 0; maxrecursion = 8;
info = 1; criterion = 1e-8; stepsize = 1; parmvary = []; minparmdiff = 1e-5;

% third argument (first varargin) could be x0, or else it could be first
% of a list of option pairs.
if nargin >= 3 && (isnumeric(varargin{1}) || isstruct(varargin{1}))
  x0 = varargin{1};
  varargin = varargin(2:end);
elseif nargin < 3 || isempty(varargin{1}) || ~isnumeric(varargin{1})  % x0 not given, get it from xstar
 x0 = get(w1, 'xstar');
 if isempty(x0)
   error('findgait requires an initial guess x0')
 end
end

% if we've been called with a numeric x0, let's convert that into
% a trial solution
if isnumeric(x0)
  solinit = bvpsetup(w1, x0);
elseif isstruct(x0) % looks like it is a structure suitable for solinit
  solinit = x0;
end

% Deal with the parameter list
options = varargin;
opt_argin = options;
while length(opt_argin) >= 2,
 opt = opt_argin{1};
 val = opt_argin{2};
 opt_argin = opt_argin(3:end);
 switch opt
   case {'criterion','NMax','solver'}
     % options for bvpsearch, do nothing
   case 'info'
     info = val;
   case 'parmvary'
     parmvary = val;
   case 'minparmdiff'
     minparmdiff = val;
   case 'maxrecursion'
     maxrecursion = val;
   case {'AbsTol','RelTol'}  
     % options for bvpsearch, do nothing
   otherwise
     error('findgaitb options: extrapolate, parmvary, info, criterion, stepsize')
 end
end

% If the user didn't specify which parameter is being varied, pick the
% parameter which undergoes the largest change
if isempty(parmvary)
 parms2 = get(w2, 'parms'); fnames = fieldnames(parms2);
 % allow for non-numeric parameters (which cannot be changed by findgaitb)
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
   error('findgaitb could not determine which parameter is being varied')
 end
 parmvary = fnames{iparm};
 options = {options{:}, 'parmvary', parmvary}; % add this knowledge to parameters
else % we know exactly which parm is being varied
 pvalues1 = get(w1, parmvary); pvalues2 = get(w2, parmvary);
%  p = abs(pvalues1-pvalues2)/max(abs(pvalues1),abs(pvalues2));
 p = abs(pvalues1-pvalues2)/max([abs(pvalues1),abs(pvalues2),any(pvalues1==0),any(pvalues2==0)]);
end

[wstar, cnvrg, sol] = findgaitsearcher(w1, w2, solinit, recursion);

if info > 0
  bcres = bvpbcfun(wstar, sol);
  fprintf('exit findgaitb  cnvrg = %d  bcres = %g\n', cnvrg, max(abs(bcres))); 
end

% Effective end of findgaitb

%% Subfunction findgaitsearcher for recursion
function [wstar, cnvrg, sol] = findgaitsearcher(w1, w2, solinit, recursion);
% start with the gait w1 and see if its fixed point works for w2

infoindent = repmat('  ', 1, recursion);

if info > 0  % diagnostic information
 fprintf(1,'start: ');
 if isempty(parmvary) % if we don't know which parm is being varied,
   display(w1); % display the whole walk object
 else
   fprintf(1, '%s = %g  ', parmvary, get(w1,parmvary));
 end
 fprintf(1,'target: ');
 if isempty(parmvary)
   display(w2);
 else
   fprintf(1, '%s = %g\n', parmvary, get(w2,parmvary));
 end
end

% The first thing to try is to converge on the target directly
try
 [wstar, cnvrg,sol] = bvpsearch(w2, solinit, options{:},'info',info-1);
catch ME
 if isempty(findstr('SingJac', ME.identifier)), % a real error
   rethrow(ME); % anything but the singular Jacobian is just an error
 end
 % singular Jacobian, meaning a failed search but not a serious error
 cnvrg = 0;
 if info > 0
   fprintf(1,'bvpsearch 1 failed\n')
 end
end

if cnvrg, % if you're already done, return successfully
 if info > 1
   fprintf('%sSucceed, findgaitb return recursion = %d\n', infoindent, recursion);
 end
 drawnow;
 return
else % not converged, get out if parameters are basically not changing
 if info >= 2, fprintf(1, 'parmdiff = %g\n',p); end
 if p < minparmdiff % there's basically no difference with the current gait
   % and we shouldn't bother with interpolating further
   fprintf(1, 'interpolation limit reached: minparmdiff\n');
   wstar = w1; 
   sol = solinit;
   return
 end   
end

if recursion >= maxrecursion
  if info >= 2,
    fprintf('max recursion limit reached in findgaitb, recursion = %d\n', recursion)
  end
  cnvrg = 0;
  sol = solinit;
  return;
end

% If direct convergence doesn't work, then it's time to interpolate
% w1 works, w2 doesn't, so let's interpolate
if info > 0, fprintf(1,'%sbisect 1 begin  ', infoindent); end
[wintermediate, cnvrg, sol] = findgaitsearcher(w1, interpolategait(w1,w2),solinit, recursion+1);
%if info > 0, fprintf(1,'%s end\n'); end

if cnvrg, % the interpolation succeeded, so try to get the target again
%           this time starting from the intermediate solution
  cnvrg = 0;
  solinit = sol;
  if info > 0, fprintf(1,'%sSucceed, bisect 2 begin  ', infoindent); end;
  [wstar, cnvrg, sol] = findgaitsearcher(wintermediate, w2, solinit, recursion+1);
else
 if info > 0, fprintf(1,'%sbisect 1 no convergence\n', infoindent); end
 wstar = wintermediate; % try to return the one that converged
end

if info > 1
  if cnvrg
    fprintf('%sSucceed cnvrg = %d\n', infoindent, cnvrg); 
  else
    fprintf('%sNo convergence cnvrg = %d\n', infoindent, cnvrg);
  end
end

end % subfunction findgaitsearcher

end