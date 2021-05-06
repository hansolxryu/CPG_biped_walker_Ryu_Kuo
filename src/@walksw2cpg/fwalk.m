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

% kappa = parms.Kp; % there's a pelvis spring

gamma = parms.gamma; Kp = parms.Kp; g = parms.g; Mp = parms.Mp; M = parms.M;
L = parms.L; C = parms.C; R = parms.R; Ip = parms.Ip; Il = parms.Il;

CutY = parms.CutY; % to simulate fictive locomotion. gradually cut. 
Lgain = parms.Lgain;

% noise = parms.testParm;
if(~isempty(parms.testParmMdl))
    noise_total = ppval(parms.testParmMdl, parms.GlobalT + t);
else
    noise_total = zeros(6,1);
end

sg = sin(gamma);
cg = cos(gamma);

%% udot_h %
q1_h = x(5); q2_h = x(6);
u1_h = x(7); u2_h = x(8);

c1 = cos(q1_h); c2 = cos(q2_h); c12 = cos(q1_h-q2_h);
s1 = sin(q1_h); s2 = sin(q2_h); s12 = sin(q1_h-q2_h);

MM_h = zeros(2,2);
rhs_h = zeros(2,1);

MM_h(1,1) = Il - 2*C*M*R + 2*C*c1*M*R - 2*L*M*R + 2*c1*L*M*R - 2*L*Mp*R + ...
    2*c1*L*Mp*R + M*(C*C) + M*(L*L) + Mp*(L*L) + 4*M*(R*R) - ...
    4*c1*M*(R*R) + 2*Mp*(R*R) - 2*c1*Mp*(R*R);

MM_h(1,2) = C*c12*L*M - C*c12*M*R + C*c2*M*R + c12*L*M*R - c2*L*M*R - ...
    c12*M*(L*L);

MM_h(2,1) = MM_h(1,2);

MM_h(2,2) = Il - 2*C*L*M + M*(C*C) + M*(L*L);

rhs_h(1) = - cg*g*M*(-C + R)*s1 - cg*g*M*(-L + R)*s1 - ...
    cg*g*Mp*(-L + R)*s1 + (-(g*M*R) + c1*g*M*(-C + R))*sg +  ...
    (-(g*M*R) + c1*g*M*(-L + R))*sg +  ...
    (-(g*Mp*R) + c1*g*Mp*(-L + R))*sg +  ...
    s1*(C*M*R + L*M*R + L*Mp*R - 2*M*(R*R) - Mp*(R*R))*  ...
    (u1_h*u1_h) + ((C*M*R - L*M*R)*s2 +  ...
    s12*(-(C*L*M) + C*M*R - L*M*R + M*(L*L)))*(u2_h*u2_h);

rhs_h(2) = s12*(C*L*M - C*M*R + L*M*R - M*(L*L))* ...
    (u1_h*u1_h) + g*(-C + L)*M*sin(gamma - q2_h);

%% udot %
q1 = x(1); q2 = x(2);
u1 = x(3); u2 = x(4);

c1 = cos(q1); c2 = cos(q2); c12 = cos(q1-q2);
s1 = sin(q1); s2 = sin(q2); s12 = sin(q1-q2);

MM = zeros(2,2); rhs = zeros(2,1);

MM(1,1) = Il - 2*C*M*R + 2*C*c1*M*R - 2*L*M*R + 2*c1*L*M*R - 2*L*Mp*R + ...
    2*c1*L*Mp*R + M*(C*C) + M*(L*L) + Mp*(L*L) + 4*M*(R*R) - ...
    4*c1*M*(R*R) + 2*Mp*(R*R) - 2*c1*Mp*(R*R);

MM(1,2) = C*c12*L*M - C*c12*M*R + C*c2*M*R + c12*L*M*R - c2*L*M*R - ...
    c12*M*(L*L);

MM(2,1) = MM(1,2);

MM(2,2) = Il - 2*C*L*M + M*(C*C) + M*(L*L); % I1 - M(C-L)^2

rhs(1) = - cg*g*M*(-C + R)*s1 - cg*g*M*(-L + R)*s1 - ...
    cg*g*Mp*(-L + R)*s1 + (-(g*M*R) + c1*g*M*(-C + R))*sg +  ...
    (-(g*M*R) + c1*g*M*(-L + R))*sg +  ...
    (-(g*Mp*R) + c1*g*Mp*(-L + R))*sg +  ...
    s1*(C*M*R + L*M*R + L*Mp*R - 2*M*(R*R) - Mp*(R*R))*  ...
    (u1*u1) + ((C*M*R - L*M*R)*s2 +  ...
    s12*(-(C*L*M) + C*M*R - L*M*R + M*(L*L)))*(u2*u2);

rhs(2) = + s12*(C*L*M - C*M*R + L*M*R - M*(L*L))* ...
    (u1*u1) + g*(-C + L)*M*sin(gamma - q2);

%%
C_matrix = [1 0 0 0; 0 1 0 0];
y = C_matrix*x(1:4) + noise_total(3:4);
y_h = C_matrix*x(5:8);

if(parms.internalState<0) % confused GC. acting on opposite leg.
    y = [y(2) y(1)]';
end

if(isinf(Lgain)) 
    % for pure FB, measurement overrides estimate.
    % torque will be generated sorely based on the measurement. 
    q1_h = y(1);
    q2_h = y(2);
end

% Torque based on estimated state %
Tst_h = parms.Tst;
Tsw_h = -(Kp*(q2_h));
Tq_h = [Tst_h; Tsw_h];
udot_h = MM_h\(rhs_h+Tq_h(1:2));

% Torque applied to the real legs, considering GC error
Tq = Tq_h;
if(parms.internalState<0) % wrong leg
    udot = MM\(rhs + flipud(Tq));
else
    udot = MM\(rhs + Tq);
end

dX = [u1; u2; udot] + [0;0;noise_total(1:2)];

if(~isinf(Lgain)) % finite L gain
    dX_h = [u1_h; u2_h; udot_h] + Lgain*(y*(1-CutY) - y_h);
    % may want to gradually cut Y to simulate fictive locomotion
    % otherwise, dX_h = dX_h dot + Lgain*(y-y_h)
else
    dX_h = dX;
    % pure FB, estimator is overriden by actual state. 
    % We need y_hat only to generate torque command.
    % for pure FB, measurement is used to generate torque command, and
    % estimator dynamics does not matter.
    % We will later figure out y_hat outside of the ode, by adding
    % noise to the actual state. 
end

if(CutY==1) % fictive locomotion
    xdot = [dX_h; dX_h];
else        % intact locomotion
    xdot =  [dX; dX_h];
end

end
