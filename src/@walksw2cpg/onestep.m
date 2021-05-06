function [xe,te,x,t,energies,match,int_match] = onestep(w, x0, varargin)
% ONESTEP   forward integration of 2-segment, simple walker for one step.
%   [xc,tc,x,t,energies] = onestep(w, x0, flag) integrates the ODE
%   until foot contact occurs, then computes momentum transfer.
%   Input is w (walksw2 object), x0 (state vector).
%   Output is xc, the state following impact; tc, the 
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
%     onestep(w, x0, 'tf', 5.5) specifies the ending time of the
%       simulation.
%   walksw2 is a simplest walker 2-D object
 
% Arthur D. Kuo, see:
% Kuo, A. D. (2002) Energetics of actively powered locomotion using the 
%   simplest walking model, Journal of Biomechanical Engineering, 124: 113-120. 
% Kuo, A. D. (2001) A simple model predicts the step length-speed 
%   relationship in human walking, Journal of Biomechanical Engineering, 
%   123: 264-269. 
 
 
anim = 0; AbsTol = 1e-7; RelTol = 1e-5; tf = 16; plotit = 0;

if nargin < 2 || (exist('x0') && isempty(x0)) % x0 not given
  % if xstar exists
  x0 = get(w, 'xstar');
  if isempty(x0)
    error('onestep requires an initial guess x0') 
  end 
end
 
ts_global = 0;
x0_h = x0; % default initial guess = x0
opt_argin = varargin;
internal = false; % internal model 
Lgain = get(w, 'Lgain');

while length(opt_argin) >= 2
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
    case 'x0_h'
      x0_h = val;
    case 'ts_global'
      ts_global = val;
    case 'internal'
      internal = val;
    otherwise
      error('onestep options: anim, AbsTol, RelTol, tf')
  end
end
 
options = odeset('events',@heelstrikeevent,'AbsTol', ...
    AbsTol,'RelTol',RelTol); %, 'MaxStep', 0.001);
 
t = [];
x = [];
state_matching = [];

ti = 0; % initial t
ts = 1/anim/64; % sampling interval

if(~anim)
    w.parms.testParm = []; 
end

terminal = false;
fall_down = false;

model_hs = isequal(Lgain, zeros(4,2));
% let estimator to have its own HS. only for pure FF.  

xe = nan;
te = nan;
energies = nan; % need to check

parms = get(w,'parms');

while(~terminal)    
    % This integrates the ode
    if(anim==0) % samping interval not specified
        yspan = [ti tf];
    else
        yspan = [ti:ts:tf];
    end
    
    internalState_thisItr = parms.internalState;
    % -1 if wrong GC; possible only for pure FF
    % will pass to fwalk

    parms.GlobalT = ts_global;    
    sol = ode45(@fwalk,yspan,[x0(:); x0_h(:)],options,w,parms);
 
    if(~isempty(sol.xe)&&sol.x(end)<0.95*tf) % some HS happened
        te = sol.xe(end)';
        xe = sol.ye(:,end)';

        if((sol.ie(end)==1)||(sol.ie(end)==2))
            if(length(sol.ie)>=2&&sol.xe(end)==sol.xe(end-1))
                % both model terminated
                last_ie = sol.ie(end-1:end); 
            else
                last_ie = sol.ie(end);
            end
            % last_ie: (1) or (2) or (1&2) or (...3)
            if(last_ie(1) == 1) % real model heelstrike
                terminal = true;     
                if(model_hs) % for pure FF
                    parms.internalState = parms.internalState * -1;
                    % will be flipped again if 2 also heel-striked.
                    % otherwise, GC state of internal model gets flipped. 
                end
            end            
            if(last_ie(end) == 2) % internal model heelstrike
                if(internal&&model_hs)
                    [xe(5:8),~] = s2stransition(xe(5:8), w);
                    % s2s transition of estimator. 
                    % applied only for pure FF. 
                end 
                x0_h = xe(5:8);
                x0 = xe(1:4);
                ti = te; % check duplication
                if(model_hs) % for pure FF
                    parms.internalState = parms.internalState * -1;
                end
            end            
        elseif(sol.ie(end)==3)
            terminal = true;
            fall_down = true;            
        end
    else % no HS happened withing given time. 
        terminal = true;
        fall_down = true;               
    end
    
    if anim  % intended for animation
        t_ = linspace(sol.x(1), sol.x(end), ...
            round((sol.x(end)-sol.x(1))*anim)+1)'; 
        % linearly interpolate desired number

        if(~isempty(t_)) 
            x_ = deval(sol, t_)';
            t = [t; t_];
            x = [x; x_];
        end
        
        state_matching = [state_matching; ...
            ones(size(t_))*internalState_thisItr];
    else
        t = [t; sol.x']; 
        x = [x; sol.y'];        
        state_matching = [state_matching; ...
            ones(size(sol.x'))*internalState_thisItr];
    end    
end
 
if(~fall_down)
    % didn't fall, terminated -> actual model had HS. 
    [xe(1:4),energies] = s2stransition(xe(1:4), w);

    if(~(internal&&model_hs)) % NOT pure FF. s2s on estimator anyways. 
        [xe(5:8),~] = s2stransition(xe(5:8), w);
    end 
end

% simulate without state estimator.
if(~internal)
    xe = xe(1:4);
    x = x(:, 1:4);
end

match = state_matching; % output state-match. 
int_match = parms.internalState;

if nargout == 0 || plotit
    plotstep(w,x,t);
end 