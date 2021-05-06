function [varargout] = energy(varargin)
% ENERGY  returns total energy of simple 2-D walker 
% E = energy(walk, x) takes in the state vector(s) and
% returns the total energy, kinetic energy, and potential
% energy of the walker for that state vector.
% E is a structure with the following variables:
%   E.total = total kinetic & potential energy
%   E.KE    = total kinetic only
%   E.PE    = total potential energy
%   E.PEg   = total gravitational potential energy
%   E.PEs   = total spring potential energy
% E = energy(walk) uses the fixed point 'xstar' for x and computes
%   energies over a full step; returns arrays of energies.
% With no outputs and x as array of state vectors, the energies
%   will be plotted.
% [KE, PE, PEg, PEs] = energy(w) returns the same information
%   as scalars, if the user does not want a structure
% Optionally can be called as energy(x, walk) for legacy code

ts = 0; % by default no time

if nargin == 1,
  walk = varargin{1};
  [xc,tc,x,ts] = onestep(walk);
elseif nargin == 2 % we are passed either walk, x or x, walk
  if isreal(varargin{1}) % looks like x, walk
    x = varargin{1};
    walk = varargin{2};
  else                   % probably walk, x
    walk = varargin{1};
    x = varargin{2};
  end
else
  error('incorrect number of arguments');
end

parms = get(walk, 'parms');
gamma = parms.gamma; Kp = parms.Kp; 

for i = 1:size(x,1) 
    q1 = x(i,1); q2 = x(i,2); u1 = x(i,3); u2 = x(i,4);

    KE(i) = 0.5*u1^2;
    PEg(i) = cos(gamma-q1);
    PEs(i) = 0; %0.5*Kp*q2^2;
    PE(i) = PEg(i) + PEs(i);
end

E.total = KE + PE;
E.KE = KE;
E.PE = PE;
E.PEg = PEg;
E.PEs = PEs;

varargout{1} = E;

if nargout == 0 && ~isvector(x) % a matrix of states; no outputs
  if isscalar(ts) % time not given
    ts = linspace(0,100,size(x,1)); % Time as a percentage of the length of the series of states
  end
  plot(ts,E.total-1,ts,KE,ts,PEg-1); legend('total','KE','PE');
elseif nargout == 4 % user wants scalar outputs instead of a structure
    varargout{1} = KE;
    varargout{2} = PE;
    varargout{3} = PEg;
    varargout{4} = PEs;
end
