function [value,isterminal,direction] = heelstrikeevent(t,x,walk,parms)
% HEELSTRIKEEVENT Returns event location info for detecting heelstrike,
% where we'll stop the simulation.
% [value, isterminal, direction] = heelstrikeevent(t, x, walk)
% where t is time, x is the state, and parms is a walk is a walk object containing
% the parameters

% Arthur D. Kuo, see:
% Kuo, A. D. (2002) Energetics of actively powered locomotion using the 
%   simplest walking model, Journal of Biomechanical Engineering, 124: 113-120. 
% Kuo, A. D. (2001) A simple model predicts the step length-speed 
%   relationship in human walking, Journal of Biomechanical Engineering, 
%   123: 264-269. 

q1 = x(1); q2 = x(2); u1 = x(3); u2 = x(4);

value = (-q2 + 2*q1); % equivalent to equal leg angles
isterminal = q1 < 0;   % stop if we've gone past midstance
direction = -1;   % and going negative
