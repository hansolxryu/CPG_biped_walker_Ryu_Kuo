function [xnew,energies] = s2stransition(x, walk, parms)
% xnew = s2stransition(x, walk) returns the new state given the state 
% prior to heel contact and a toe-off impulse P.  The effect of the
% toe-off impulse and the heel impact are computed as perfectly inelastic 
% impulses, and impulse-momentum and conservation of angular momentum
% principles are used to find the velocities after impact. There is no 
% change in configuration except that after heel strike the legs are 
% switched.  For heel strike, a component of angular momentum conserved 
% about the point of contact, and for the trailing leg about the hip.
% walk is a walksw2 object containing parameters.
% [xnew,energies] = s2stransition(x, walk) also returns the energies
% before and after the step-to-step transition
% x = [qstance, qsswing ustance uswing]'
% where qsswing is defined relative to the stance leg

% Arthur D. Kuo, see:
% Kuo, A. D. (2002) Energetics of actively powered locomotion using the 
%   simplest walking model, Journal of Biomechanical Engineering, 124: 113-120. 
% Kuo, A. D. (2001) A simple model predicts the step length-speed 
%   relationship in human walking, Journal of Biomechanical Engineering, 
%   123: 264-269. 

% CHANGES
%   Added parms to input list to speed up simulations

if nargin < 3
  parms = get(walk, 'parms');
end
P = parms.P; % the toe-off impulse

q1 = x(1); q2 = x(2); u1 = x(3); u2 = x(4);

c2t = cos(2*q1); s2t = sin(2*q1);

xnew(1,1) = -q1;
xnew(1,2) = -2*q1;
xnew(1,3) = c2t*u1 + s2t*P;
xnew(1,4) = c2t*(1-c2t)*u1 + (1-c2t)*s2t*P;

%mvcombefore = [-u1*cos(q1) -u1*sin(q1)];
%mvcomintermediate = mvcombefore + P*[-sin(q1) cos(q1)];
%mvcomafter = [-xnew(3)*cos(q1) -xnew(3)*sin(q1)];

energies.before = energy(walk, x); % energies before impulse
energies.intermediate = energies.before; energies.intermediate.KE = energies.intermediate.KE + 0.5*P*P; 
energies.intermediate.total = energies.intermediate.total + 0.5*P*P;
energies.after = energy(walk, xnew);
energies.pushoffwork = 0.5*P*P;
energies.heelstrikework = energies.after.KE - energies.intermediate.KE;
