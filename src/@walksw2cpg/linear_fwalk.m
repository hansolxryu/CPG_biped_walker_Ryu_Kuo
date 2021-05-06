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
A0 = parms.A0;
A = parms.A0;

Lgain = parms.Lgain;

% noise = parms.testParm;
if(~isempty(parms.testParmMdl)&&~isfield(parms, 'discreteU'))
%     noise_t = noise(:,1);
%     noise_total = interp1(noise(:,1),noise(:,2:end),...
%         parms.GlobalT + t,'spline')';
    noise_total = ppval(parms.testParmMdl, parms.GlobalT + t);
%     noise_total = zeros(6,1)+randn(6,1)*0.00001;
else
    noise_total = zeros(6,1);
end

C_matrix = [1 0 0 0; 0 1 0 0];
y = C_matrix*x(1:4) + noise_total(3:4);
y_h = C_matrix*x(5:8);


dX_h = A*x(5:8) + Lgain*(y - y_h);
dX = A*x(1:4) + [0;0;noise_total(1:2)];

xdot =  [dX+[0;0;0.3*sin(t); 0.5*sin(t)];...
    dX_h+[0;0;0.3*sin(t); 0.5*sin(t)]];
end
