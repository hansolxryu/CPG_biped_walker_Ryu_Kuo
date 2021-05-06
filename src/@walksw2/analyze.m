function [grf, grt, cop, pcom] = analyze(walk,x0)
% analyze   Computes ground reaction forces and COM work for walking
%           model.  
% [grf, grt, cop, pcom] = analyze(t, x, walk) computes the ground
% reaction forces in grf, ground reaction torques in grt, and com work rate in pcom, given input
% consisting of time t, state x, and a walk object containing the
% parameters.
% For this model, grf = [fx fy], the forward horizontal and vertical 
% components.  grt = tz, the moment about the ground3 axis.
% cop = copx, the center of pressure along the forward horizontal axis.

if nargin == 0
  error('analyze: need a walk object as first argument');
elseif nargin == 1 % optional arguments
  if ~isa(walk, 'walksw2') % assumes walksw2 or child thereof
        error('analyze: need a walk object as first argument');
  end
  [xe,te,x,t] = onestep(walk);
elseif nargin == 2
  [xe,te,x,t] = onestep(walk,x0);
elseif nargin > 2
  error('incorrect number of arguments');
end

[grf,grt,cop] = groundreactionforces(walk, x, t);
%[KE,PE,PEg,PEs] = energy(walk, x);
e = energy(walk, x);

if nargout == 0
    clf
    subplot(2,1,1);
    plot(t,grf(:,1),t,grf(:,2)); legend('Fx','Fy');
    subplot(2,1,2);
    plot(t,e.KE,t,e.PE-1,t,e.PEg-1,t,e.PEs); legend('KE','PE','PEg','PEs');
end	

pcom = []; % No center-of-mass work rate in simplest walking model