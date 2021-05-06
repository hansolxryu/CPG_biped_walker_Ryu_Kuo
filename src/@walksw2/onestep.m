function [xnext,te,x,t,energies] = onestep(w, x0, varargin)
% ONESTEP   simulation of 2-segment, simplest walking model for one step.
%   [xnext,tc,x,t,energies] = onestep(w, x0, flag) integrates the ODE
%   until foot contact occurs, then computes momentum transfer.
%   Input is w (walksw2 object), x0 (state vector).
%   Output is xnext, the state following impact (initial condition
%   for the next step); tc, the 
%   time of contact; t, the time vector; x, the state
%   vector over time of the simulation; energies contain the
%   energies before and after push-off and heelstrike.
%   The initial condition x0 = [qstance qswing ustance uswing]
%   and similarly, xc, both being row vectors.
%   Options: 
%     onestep(w, x0, 'anim', Nframes) returns the states
%       evenly spaced in time, with a total of Nframes per step
%     onestep(w, x0, 'AbsTol', 1e-9) specifies the absolute tolerance
%       for the integration, similarly 'RelTol' relative tolerance
%     onestep(w, x0, 'tf', 5.5) specifies the max ending time of the
%       simulation (if no event detected).
%   w is a walksw2 simplest walking model 2-D object

% Kuo, A. D. (2002) Energetics of actively powered locomotion using the 
%   simplest walking model, Journal of Biomechanical Engineering, 124: 113-120. 
% Kuo, A. D. (2001) A simple model predicts the step length-speed 
%   relationship in human walking, Journal of Biomechanical Engineering, 
%   123: 264-269. 

anim = 0; AbsTol = 1e-7; RelTol = 1e-5; tf = 5.5; plotit = 0;

if nargin < 2 || (exist('x0') && isempty(x0)) % x0 not given
  % if xstar exists
  x0 = get(w, 'xstar');
  if isempty(x0)
    error('onestep requires an initial guess x0') 
  end 
end

opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'anim'
      anim = val;
    case 'AbsTol'
      AbsTol = val;
    case 'RelTol'
      RelTol = val;
    case 'tf'
      tf = val;
    case 'plotstep'
      plotit = val;
    otherwise
      error('onestep options: anim, AbsTol, RelTol, tf')
  end
end

% detect any heelstrike past mid-stance (to avoid detecting scuffing)
options = odeset('events',@heelstrikeevent,'AbsTol',AbsTol,'RelTol',RelTol);

% This integrates the ode
sol = ode45(@fwalk,[0 tf],x0(:),options,w,get(w,'parms'));
t = sol.x'; x = sol.y';

if anim  % intended for animation
  t = linspace(0, sol.xe(end), anim+1)'; % linearly interpolate desired number
  t(end) = [];                           % of frames, plus one
  x = deval(sol, t)';
end

te = sol.xe(end)'; 
xe = sol.ye(:,end)'; 

% Now I have to do the switch:
[xnext,energies] = s2stransition(xe, w);

if nargout == 0 || plotit % plot angles if no output, or if asked for
    plotstep(w,x,t);
end	

