function [grf,grt,cop] = groundreactionforces(walk, x, t)
% GROUNDREACTIONFORCES   Computes ground reaction forces for simple
% 2-D passive walking.  
% [grf, grt, cop] = groundreactionforces(walk, x, t) computes the ground
% reaction forces in grf, and ground reaction torques in grt, given input
% consisting of time t, state x, and a walk object containing the
% parameters.
% x may be an nx4 matrix of states with nx1 vector of times t 
% or a single state vector at time t
% For this model, grf = [fx fy], the forward horizontal and vertical 
% components.  grt = tz, the moment about the ground3 axis.
% cop = copx, the center of pressure along the forward horizontal axis.

if nargin == 0
  error('groundreactionforces: need a walk object as first argument');
elseif nargin == 1 % optional arguments
  if ~isa(walk, 'walksw2')
        error('groundreactionforces: need a walk object as first argument');
  end
  [xe,te,x,t] = onestep(walk);
elseif nargin == 2
  if length(x) == 4   % and it's a vector, meaning an initial condition
    t = 0;
  elseif size(x,1) > 1 && size(x,2) == 4 % it's a matrix of states
    t = linspace(0,100,size(x,1)); % Time as a percentage of the length of the series of states
  else
    error('groundreactionforces: unknown second argument')
  end
elseif nargin > 3
  error('incorrect number of arguments');
end

parms = get(walk,'parms');

for i = 1:length(t)
    q1 = x(i,1); q2 = x(i,2); u1 = x(i,3); u2 = x(i,4);
    c1 = cos(q1); s1 = sin(q1);

    xdot = fwalk(t(i), x(i,:), walk, parms); % need state-derivative
    u1dot = xdot(3); u2dot = xdot(4);

    grf(i,1:2) = [ -c1*u1dot+s1*u1*u1 -s1*u1dot-c1*u1*u1+1 ];

    cop(i) = 0; 

    grt(i) = 0; % ground reaction torque about z axis
end

if nargout == 0 && length(t) > 1
    plot(t,grf(:,1),t,grf(:,2)); legend('Fx','Fy');
end