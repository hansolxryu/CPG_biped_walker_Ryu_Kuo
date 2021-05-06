function xdot = fwalk(t,x,walk,parms)
% FWALK   Right-hand side for Matlab ode45, simple walker 2-D
% xdot = fwalk(t, x, walk) computes the right-hand-side of the
% differential equations of motion, given time t and state x.
% walk is a walk object containing the parameters

% Arthur D. Kuo, see:
% Kuo, A. D. (2002) Energetics of actively powered locomotion using the 
%   simplest walking model, Journal of Biomechanical Engineering, 124: 113-120. 
% Kuo, A. D. (2001) A simple model predicts the step length-speed 
%   relationship in human walking, Journal of Biomechanical Engineering, 
%   123: 264-269. 

% CHANGES
%   Vectorized by Art (6/28/09) to make some codes run much faster

kappa = parms.Kp; % there's a pelvis spring
gamma = parms.gamma;

cg = cos(gamma); sg = sin(gamma);

if all(size(x) > 1) % vectorized states are in columns
  q1 = x(1,:); q2 = x(2,:); u1 = x(3,:); u2 = x(4,:);
else % just one column vector of the state
  q1 = x(1); q2 = x(2); u1 = x(3); u2 = x(4);
end

% where q2 is the swing angle relative to stance, to find
% the swing angle relative to vertical, use q1-q2

u1dot = sin(q1-gamma);  % + T
u2dot = u1dot + (u1.*u1 - cos(q1-gamma)).*sin(q2) - kappa*q2;
q1dot = u1;
q2dot = u2;
xdot = [q1dot; q2dot; u1dot; u2dot];
