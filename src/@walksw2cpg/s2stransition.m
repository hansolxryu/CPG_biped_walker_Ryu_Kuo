function [xnew,energies] = s2stransition(x, walk, parms)
% [xnew, energies] = s2stransition(x, walk) returns the new state given the state 
% prior to heel contact and a toe-off impulse P.  The effect of the
% toe-off impulse and the heel impact are computed as perfectly inelastic 
% impulses, and impulse-momentum and conservation of angular momentum
% principles are used to find the velocities after impact. There is no 
% change in configuration except that after heel strike the legs are 
% switched.  For heel strike, a component of angular momentum conserved 
% about the point of contact, and for the trailing leg about the hip.
% walk is a walk2 object containing parameters.
% x = [qstance, qsswing ustance uswing]'

% Arthur D. Kuo, see:
% Kuo, A. D. (2002) Energetics of actively powered locomotion using the 
%   simplest walking model, Journal of Biomechanical Engineering, 124: 113-120. 
% Kuo, A. D. (2001) A simple model predicts the step length-speed 
%   relationship in human walking, Journal of Biomechanical Engineering, 
%   123: 264-269. 
% McGeer (1990) Passive dynamic walking, Intl. J. Robotics Research,
%   9: 62-82.

% CHANGES
%   Added parms to input list to speed up simulations (Art, 6/2009)

if nargin < 3 % didn't receive parms in input list
  parms = get(walk, 'parms');
end
gamma = parms.gamma; Kp = parms.Kp; g = parms.g; Mp = parms.Mp; M = parms.M;
L = parms.L; C = parms.C; R = parms.R; Ip = parms.Ip; Il = parms.Il;
P = parms.P;

% sg = sin(gamma); cg = cos(gamma);

q1 = x(1); q2 = x(2); u1 = x(3); u2 = x(4);
 
c1 = cos(q1); c2 = cos(q2); c12 = cos(q1-q2);
s1 = sin(q1); s2 = sin(q2); s12 = sin(q1-q2);
c12h = cos((q1-q2)/2);

angle = atan2(-(L-R)*sin(q1), R+(L-R)*cos(q1));
Ix = P*sin(angle); Iy = P*cos(angle);


% q1 = x(1); q2 = x(2); u1 = x(3); u2 = x(4);
% c1 = cos(q1); c2 = cos(q2); c12 = cos(q1-q2);
% s1 = sin(q1); s2 = sin(q2); s12 = sin(q1-q2);

% All of this is the intermediate state of the system
% angular momentum of entire system about com, in terms of old u's
% H(t0) = H(t-) + r1 x impulse using u(t-) (just like above)
amsciold = Ix*(R+(C*M+L*(M+Mp)-R)*c1 + (C-L)*M*c2) + ...
  Iy*((C*M+L*(M+Mp)-R)*s1 + (C-L)*M*s2) + ...
  (Il+(C-L)^2*M*(M+Mp)-(C-L)^2*M*M*c12)*u1 + ...
  (Il-(C-L)^2*M*(-M-Mp+M*c12))*u2;
% angular momentum of trailing leg about hip, in terms of old u's
% Htl(t0) = Htl(t-) + r x impulse (r vector from contact pt to hip)
amthiold = -Ix*(-R+(-L+R)*c1) + Iy*(L-R)*s1 + (Il+(C-L)*M*(C-R)+(C-L)*M*R*c1)*u1;
% angular momentum of leading leg about hip, in terms of old u's
% Hll(t0) = Hll(t-)
amlhiold = (C-L)*M*((L-R)*c12 + R*c2)*u1 + (Il+(C-L)^2*M)*u2;
% linear momentum, in terms of old u's
% L(t0) = L(t-) + impulse
lmioldx = Ix - R*u1 - (C*M+L*(M+Mp)-R)*c1*u1 - (C-L)*M*c2*u2;
lmioldy = Iy - (C*M+L*(M+Mp)-R)*s1*u1 - (C-L)*M*s2*u2;
KEstab = 1/2*M*((C^2-2*C*R+2*R^2+2*(C-R)*R*c1)*u1^2) + 1/2*Il*u1^2;
KEswib = 1/2*M*((L^2-2*L*R+2*R^2+2*(L-R)*R*c1)*u1^2 + ...
  2*u1*((C-L)*((L-R)*c12+R*c2)*u2) + (C-L)^2*u2^2) + 1/2*Il*u2^2;
KEpelvb = 1/2*Mp*((L^2-2*L*R+2*R^2+2*(L-R)*R*c1)*u1^2);

% Intermediate u's, old q's 
% Row 1: amsciint (angular mom of whole system about com using interm. u's
MM(1,1) = (Il+(C-L)^2*M*(M+Mp)-(C-L)^2*M*M*c12);
MM(1,2) = MM(1,1);
MM(1,3) = 0; MM(1,4) = 0;
% Row 2: amthiint (ang mom of trailing leg about hip using interm. u's)
MM(2,1) = (Il+(C-L)*M*(C-R) + (C-L)*M*R*c1); MM(2,2) = 0;
MM(2,3) = (-C+L)*M*c1; MM(2,4) = -(C-L)*M*s1;
% Row 3: lmiintx (linear mom x component, interm. u's)
MM(3,1) = (-R+(-C*M-L*(M+Mp)+R)*c1);
MM(3,2) = (-C+L)*M*c2; MM(3,3) = 1; MM(3,4) = 0;
% Row 4: lmiinty (linear mom y component, interm. u's)
MM(4,1) = (-C*M-L*(M+Mp)+R)*s1; MM(4,2) = (-C+L)*M*s2;
MM(4,3) = 0; MM(4,4) = 1;
% Solve for full intermediate state of system, where u3 and u4 refer to
% the velocity of the foot: u[3] ground[1] + u[4] ground[2]
% MM * uintermediate = [amsciold; amthiold; lmioldx; lmioldy]
uintr = MM \ [amsciold; amthiold; lmioldx; lmioldy];
u1 = uintr(1); u2 = uintr(2); u3 = uintr(3); u4 = uintr(4);
amlhiint = (C-L)*M*((L-R)*c12 + R*c2)*u1 + ...
  (Il+(C-L)^2*M)*u2 + (-C+L)*M*c2*u3 - (C-L)*M*s2*u4;
% now we can compute the kinetic energy in the intermediate state
KEstai = 1/2*M*((C^2-2*C*R+2*R^2+2*(C-R)*R*c1)*u1^2 + u3^2 + u4^2 + ...
  2*u1*((-R+(-C+R)*c1)*u3 + (-C+R)*s1*u4)) + 1/2*Il*u1^2;
KEswii = 1/2*M*((L^2-2*L*R+2*R^2+2*(L-R)*R*c1)*u1^2 + (C-L)^2*u2^2 + u3^2 + u4^2 + ...
  2*u1*((C-L)*((L-R)*c12+R*c2)*u2 + (-R + (-L+R)*c1)*u3 + ...
  (-L+R)*s1*u4) - 2*(C-L)*u2*(c2*u3+s2*u4)) + 1/2*Il*u2^2;
KEpelvi = 1/2*Mp*((L^2-2*L*R+2*R^2+2*(L-R)*R*c1)*u1^2 + u3^2 + u4^2 + ...
  2*u1*((-R + (-L+R)*c1)*u3 + (-L+R)*s1*u4));
e1 = KEstai + KEswii + KEpelvi;
% u3 and u4 contain velocity of trailing foot after impulse,
% and ux and uy contain velocity of leading foot after impulse
ux = (-R+(-L+R)*c1)*u1 + (R+(L-R)*c2)*u2 + u3; % velocity of trailing foot
uy = -(L-R)*s1*u1 + (L-R)*s2*u2 + u4;          % after impulse

% New u's, new q's (post-heelstrike info)
c1 = cos(q2); c2 = cos(q1); s1 = sin(q2); s2 = sin(q1);
% The MM matrix below is used to solve for u1(t+), u2(t+), and cx and cy
% the impact collision vector
MM(1,1) = (Il+(C-L)^2*M*(M+Mp)-(C-L)^2*M*M*c12); % ang mom of whole system about com
MM(1,2) = MM(1,1);
MM(1,3) = -(R+(C-L)*M*c2 + (C*M+L*(M+Mp)-R)*c1); % MM*[u1;u2;cx;cy]
MM(1,4) = -((C-L)*M*s2 + (C*M+L*(M+Mp)-R)*s1);
MM(2,1) = (C-L)*M*((L-R)*c12 + R*c2);
MM(2,2) = Il+(C-L)^2*M;
MM(2,3) = 0; MM(2,4) = 0;
MM(3,1) = -(R+(C*M+L*(M+Mp)-R)*c1);
MM(3,2) = (-C+L)*M*c2;
MM(3,3) = -1; MM(3,4) = 0;
MM(4,1) = (-C*M-L*(M+Mp)+R)*s1;
MM(4,2) = (-C+L)*M*s2;
MM(4,3) = 0; MM(4,4) = -1;

% solve for new u's, plus collision vector cx and cy
unew = MM \ [amsciold;amthiold;lmioldx;lmioldy]; 
u1 = unew(1); u2 = unew(2); cx = unew(3); cy = unew(4);
% angular momentum of entire system about com, in terms of new u's
amscinew = cx*(-R-(C*M+L*(M+Mp)-R)*c1 - (C-L)*M*c2) + ...
  cy*(-(C*M+L*(M+Mp)-R)*s1 + (-C+L)*M*s2) + ...
  (Il+(C-L)^2*M*(M+Mp)-(C-L)^2*M*M*c12)*u1 + ...
  (Il-(C-L)^2*M*(-M-Mp+M*c12))*u2;
% angular momentum of trailing leg about hip, in terms of new u's
amthinew = (C-L)*M*((L-R)*c12 + R*c2)*u1 + (Il+(C-L)^2*M)*u2;
% angular momentum of leading leg about hip, in terms of new u's
amlhinew = cx*(-R+(-L+R)*c1) - cy*(L-R)*s1 + (Il+(C-L)*M*(C-R)+(C-L)*M*R*c1)*u1;
% linear momentum, in terms of new u's
lminewx = -cx - R*u1 - (C*M+L*(M+Mp)-R)*c1*u1 - (C-L)*M*c2*u2;
lminewy = -cy - (C*M+L*(M+Mp)-R)*s1*u1 - (C-L)*M*s2*u2;

KEstaa = 1/2*M*((C^2-2*C*R+2*R^2+2*(C-R)*R*c1)*u1^2) + 1/2*Il*u1^2;
KEswia = 1/2*M*((L^2-2*L*R+2*R^2+2*(L-R)*R*c1)*u1^2 + ...
  2*u1*((C-L)*((L-R)*c12+R*c2)*u2) + (C-L)^2*u2^2) + 1/2*Il*u2^2;
KEpelva = 1/2*Mp*((L^2-2*L*R+2*R^2+2*(L-R)*R*c1)*u1^2);

KEb = KEstab+KEswib+KEpelvb;
KEi = KEstai+KEswii+KEpelvi;
KEa = KEstaa+KEswia+KEpelva;

% 
xnew(1) = x(2);
xnew(2) = x(1);
xnew(3) = unew(1);
xnew(4) = unew(2);

energies.before = energy(x, walk); % energies before impulse
energies.intermediate = energies.before; energies.intermediate.KE = KEi; energies.intermediate.total = KEi + energies.intermediate.PE;
energies.after = energy(xnew, walk);
energies.pushoffwork = 0.5*(Ix*uintr(3)+Iy*uintr(4));
energies.heelstrikework = 0.5*(cx*ux+cy*uy);

energies.impulse.pushoff = [Ix; Iy];
energies.impulse.collision = [cx; cy];
