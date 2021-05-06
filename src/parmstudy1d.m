function [ws,parmrangesucc,cnvrg] = parmstudy1d(w, parmrange, parmname, varargin)
% ws = parmstudy1d(w, parmrange, 'parmname') performs a one-dimensional
%   parameter study of the class w, setting the parameter named in 'parmname'
%   to the values in the array parmrange. The output argument ws is an array
%   of walking objects, which are referred to be index, e.g. ws(1).
% [ws, parmrangesucc, cnvrg] = parmstudy1d... also returns two arrays, 
%   parmrangesucc being the parameter range that was actually successful, and 
%   cnvrg an array of same size, with 1 meaning
%   successful convergence and 0 meaning a failure. In cases of "failure"
%   there should be a gait at the extreme of the parameter range, just not
%   at the intended value. So cnvrg(1) and cnvrg(end) should indicate
%   whether either end of desired range was not met, and parmrangesucc(1)
%   and end should indicate what those parameters were.  
% Note that 

% Added by Art July 2008

% CHANGES
%   Modified by Art to include an array of cnvrg indicators of success
%   (6/2009)
%   Modified by Shawn to include parmrangesucc (the succ suffix added
%   by Art

% Example: This performs a parameter study varying gravity 'g' from
% 0.16 - 2, and returns an array of walking objects
%   ws = parmstudy1d(walk2('normal'), linspace(0.16, 2, 20), 'g');

p0 = get(w, parmname); x0 = get(w, 'xstar');

% We start with a w that is presumably somewhere within parmrange. 
% Use lowindices and highindices to break the parms into two ranges, 
% below and above the starting gait.
lowindices = find(p0 > parmrange);
highindices = find(p0 < parmrange);
% put the known w smack in the middle
myparmrange = [parmrange(lowindices) p0 parmrange(highindices)];
index0 = length(lowindices) + 1;
ws(index0) = w; % put the starting point in the array
cnvrg = zeros(size(parmrange)); cnvrg(index0) = 1;
ibegin = index0; iend = index0; % try to keep track of the gaits that fail completely

% Descending ladder
for i = fliplr(lowindices)
  try 
    [ws(i),cnvrg(i)] = findgait(w, set(w, parmname, parmrange(i)), x0, 'parmvary', parmname, 'extrapolate', 1, varargin{:});
  catch
    ibegin = i+1;
    break
  end
  w = ws(i); 
  myparmrange(i) = get(w, parmname); % store the actual value that was used
  if ~cnvrg, break; end;  % give up if we didn't converge
  ibegin = i;
  x0 = get(w, 'xstar');
  if i > 1 % use whatever succeeded to try to extrapolate the next gait
    x0 = extrapolategait(ws(i:min(i+2,index0)), set(w, parmname, parmrange(i-1)), 'parmvary', parmname);
  end
end

% After finishing the descending ladder, 
% try to extrapolate ahead for the ascending ladder
if index0 > 1 && length(highindices) > 1
  x0 = extrapolategait(ws(max(1,index0-2):index0), set(w, parmname, parmrange(index0+1)), 'parmvary', parmname);
end

% Ascending ladder
for i = highindices
  try
    [ws(i),cnvrg(i)] = findgait(w, set(w, parmname, parmrange(i)), x0, 'parmvary', parmname, 'extrapolate', 1, varargin{:});
  catch
    iend = i-1;
    break
  end
  w = ws(i); 
  myparmrange(i) = get(w, parmname); % store the actual value that was used
  if ~cnvrg, break; end;  % give up if we didn't converge
  iend = i;
  x0 = get(w, 'xstar');
  if i < length(parmrange)
    x0 = extrapolategait(ws(max(1,i-2):i), set(w, parmname, parmrange(i+1)), 'parmvary', parmname);
  end
end

ws = ws(ibegin:iend); % return only usable gaits, although not necessarily
% the intended parameter values
parmrangesucc = myparmrange(ibegin:iend); % the actual values used