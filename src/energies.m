function [E,xs,ts,deltaE] = energies(varargin)
% ENERGIES returns and/or plots the energies for an entire step
%   e = energies(w) performs a simulation of walk object w and returns the
%     energies over time in the structure of arrays e.
%   e = energies(xs, ts, w) uses the time and state vectors provided.
%   energies(xs, ts, w) without an output argument automatically 
%     yields a plot of all of the energy components together
%   [e, xs, ts] = ... also returns the state and time vectors
%
% energies calls the energy function, which usually returns a structure
% such as e.total, e.KE, e.PE, etc. The present function returns the same
% fields, but rearranged into a structure of arrays to simplify plotting.
% For example, it will return e.total which is an array of energies over
% the time of a whole step.
%
% Extra output [E, xs, ts, deltaE] = ... includes the max energy change
% in the first field of E (usually the total energy).

% Art Kuo 7/2012

if nargin == 1 % then all we have is a walk object
  w = varargin{1};
  [xe, te, xs, ts] = onestep(w);
elseif nargin == 2 % then we have xs and w, but not time
  xs = varargin{1};
  w = varargin{2};
  ts = 1:size(xs,1); % just make time a sample count
elseif nargin == 3 % then we have ts, xs, and w
  xs = varargin{1};
  ts = varargin{2};
  w = varargin{3};
end

N = size(xs, 1);
estruct = energy(xs(1,:), w); % get one sample of the energy structure
names = fieldnames(estruct);
numNames = length(names);
energyarray = zeros(N, numNames); % make a single array to store all the
  % energy components
  
for i = 1:N % loop over time steps
  estruct = energy(xs(i,:), w); % get one sample of the energy structure
  for j = 1:numNames % and get each field in the structure
    energyarray(i,j) = estruct.(names{j});
  end
end

% reconstruct everything into a structure of arrays with appropriate
% fieldnames, only if there is an output argument
if nargout > 0
  for j = 1:numNames
    E.(names{j}) = energyarray(:,j);
  end
  deltaE = max(energyarray(:,1))-min(energyarray(:,1));
else
  plot(ts, energyarray)
  xlabel('time'); ylabel('energy');
  legend(names)
end

