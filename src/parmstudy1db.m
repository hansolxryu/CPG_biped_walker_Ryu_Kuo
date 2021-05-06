function [ws,cnvrg, sols] = parmstudy1d(w, parmrange, parmname, varargin);
% ws = parmstudy1db(w, parmrange, 'parmname') performs a one-dimensional
%   parameter study of the class w, setting the parameter named in 'parmname'
%   to the values in the array parmrange. The output argument ws is an array
%   of walking objects, which are referred to be index, e.g. ws(1).
% [ws, cnvrgs, sols] = parmstudy1db... returns array cnvrg with 1's meaning
%   a particular gait was found successfully.
%
% This version uses boundary value problem solver findgaitb

% Added by Art Jun 2009

% Example: This performs a parameter study varying gravity 'g' from
% 0.16 - 2, and returns an array of walking objects
%   ws = parmstudy1d(walk2('normal'), linspace(0.16, 2, 20), 'g');

p0 = get(w, parmname); x0 = get(w, 'xstar');
solinit1 = bvpsetup(w);

lowindices = find(p0 > parmrange);
highindices = find(p0 < parmrange);
% put the known w smack in the middle
myparmrange = [parmrange(lowindices) p0 parmrange(highindices)];
index0 = length(lowindices) + 1;
ws(index0) = w; % put the starting point in the array
cnvrg = zeros(size(parmrange)); cnvrg(index0) = 1;

solinit = solinit1;

% first let's descend
for i = fliplr(lowindices)
  try 
    [ws(i),cnvrg(i),errmsg,wsuccess,solinit] = ...
      findgaitb(w, set(w, parmname, parmrange(i)), solinit, 'parmvary', parmname, varargin{:});
  catch ME
    warning('parmstudy1b failed on low search');
    rethrow(ME)
  end
  if ~cnvrg(i), warning('could not converge in parmstudy1db'); break; end
  w = ws(i); x0 = get(w, 'xstar');
  if i > 1 
    %x0 = extrapolategait(ws(i:min(i+2,index0)), set(w, parmname, parmrange(i-1)), 'parmvary', parmname);
  end
end

% and try to extrapolate ahead
if index0 > 1 && length(highindices) > 1
  %x0 = extrapolategait(ws(max(1,index0-2):index0), set(w, parmname, parmrange(index0+1)), 'parmvary', parmname);
end

solinit = solinit1;

for i = highindices
  try
    [ws(i),cnvrg(i),errmsg,wsuccess,solinit] = ...
      findgaitb(w, set(w, parmname, parmrange(i)), solinit, 'parmvary', parmname, varargin{:});
  catch ME
    warning('parmstudy1b failed on high search');
    rethrow(ME)
  end
  if ~cnvrg(i), warning('could not converge in parmstudy1db'); break; end
  w = ws(i); x0 = get(w, 'xstar');
  if i < length(parmrange)
    %x0 = extrapolategait(ws(max(1,i-2):i), set(w, parmname, parmrange(i+1)), 'parmvary', parmname);
  end
end
  