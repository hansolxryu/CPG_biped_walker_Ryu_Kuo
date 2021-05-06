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

% CHANGES
%   Modified by Art to remove the ugly kludge by which we biased things for
%   long-period gaits, and avoidance of mid-stance scuffing.

% we don't need any parms, so this is currently unused
%if nargin < 4 % we are not given parameters
  %parms = get(walk, 'parms');
%end
  
% DEPRECATED:
%q10 = (x0(1) + parms.gamma); % this is an ugly kludge to avoid stumbling
% if the swing foot scuffs the ground. It is safer in future to make the
% modification to isterminal, and leave value alone.

q1 = x(1); q2 = x(2); u1 = x(3); u2 = x(4);
q1_h = x(5); q2_h = x(6); u1_h = x(7); u2_h = x(8);

if nargin < 4 % didn't receive parms, so get it from walk
    parms = get(walk, 'parms');
end
% DEPRECATED:
%x0 = get(walk, 'xstar');
%q10 = x0(1); u10 = x0(3);
L = parms.L; R = parms.R;

value = [ q1+q2; %((cos(q1)*(L-R) + cos(q2)*(-L+R)));
    q1_h+q2_h;
    (pi/2-q1)*(q1+pi/2);
    (pi/2-q1_h)*(q1_h+pi/2)*0+1];

x0 = get(walk, 'xstar');
CutY = get(walk, 'CutY');

isterminal = [(q1<-0.1*x0(1))||(CutY==1&&q1<0); 
    (q1_h<-0.1*x0(1))||(CutY==1&&q1<0); 
    true; 
    false];

direction = [-1; -1; -1; -1];  % and going negative
end