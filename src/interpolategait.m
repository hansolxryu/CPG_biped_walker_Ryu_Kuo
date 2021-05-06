function wnew = interpolategait(w1, w2, c)
% wnew = interpolategait(w1, w2) interpolates between two gaits, 
% effectively returning a gait like wnew = c*w1 + (1-c)*w2, 
% based on all parameters of interest.

% Changes
%  v 1.1  Added non-numeric parameter handling (July 6 2008)

if nargin < 3
  c = 0.5;
end

wnew = w1;

parms = get(w1,'parms');
names = fieldnames(parms);
for i = 1:length(names)
  nameparms{2*i-1} = names{i};
  if isnumeric(get(w1,(names{i})));
      nameparms{2*i} = c*get(w1,(names{i})) + (1-c)*get(w2,(names{i}));
  else
      nameparms{2*i} = get(w2,(names{i}));
  end
 
%  name = names{i};
%  parms.(name) = c*get(w1,(name)) + (1-c)*get(w2,(name));
end

%wnew = set(wnew, 'parms', parms);
wnew = set(wnew, nameparms{:});
%%%%