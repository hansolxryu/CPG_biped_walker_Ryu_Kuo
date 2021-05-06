function [xd,td,xs,ts,stepindices] = manystep(w, x0, varargin)
% MANYSTEP   integrates walking model over numsteps
%     manystep(w, x0) integrates the walking model described by w.
%   Outputs: [xd, td, xs, ts, stepindices] where xd, td are the states
%    and times at end of each step, and xs, ts are the trajectories
%    over time. Stepindices is useful as an input to animate, showing
%    which entries in xs, ts delineate steps.
%   Options: 
%     manystep(w, x0, 'numsteps', N) automatically repeats the simulation N times
%       (default 3).
%     manystep(w, x0, 'anim', Nframes) returns the states
%       evenly spaced in time, with a total of Nframes per step
%     manystep(w, x0, 'AbsTol', 1e-9) specifies the absolute tolerance
%       for the integration, similarly 'RelTol' relative tolerance
%     manystep(w, x0, 'tf', 5.5) specifies the ending time of the
%       simulation for each step.
%     manystep(w, x0, 'plotstep', 1) adds a plot in the current figure
%       (default behavior is to plot when there are no output arguments
%       but otherwise not)

% CHANGES
%   Art added 'plotstep' option 6/2009

numsteps = 3; anim = 0; doplot = 0; info = 0;

if nargin < 2 || (exist('x0','var') && isempty(x0)) % x0 not given
  % if xstar exists
  x0 = get(w, 'xstar');
  if isempty(x0)
    error('onestep requires an initial guess x0') 
  end 
end

% Deal with the parameter list
onestepoptions = {};
opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'anim'
      anim = val;
      onestepoptions = {onestepoptions{:}, opt, val}; 
    case 'fps'
      onestepoptions = {onestepoptions{:}, opt, val};
    case {'AbsTol', 'RelTol', 'tf'}
      onestepoptions = {onestepoptions{:}, opt, val};       
    case 'numsteps'
      numsteps = val;
    case 'plotstep'
      doplot = val;
    case 'info' % set this to 1 if you want to see steps counted out
      info = val; 
    otherwise
      error('onestep options: anim, AbsTol, RelTol, tf, numsteps')
  end
end

xs = []; ts = []; xdtemp = []; td = []; stepindices = []; ti = 0; 

for i = 1:numsteps
    if info, fprintf('%d ', i); end;
    [xe,te,x,t] = onestep(w, x0, onestepoptions{:});
    xs = [xs;x];  ts = [ts;ti+t]; 
    xdtemp = [xdtemp;xe]; td = [td;ti+te];
    x0 = xe;
    stepindices = [stepindices length(t)];
    ti = ti+te;
end
if info, fprintf('\n'); end

if anim    
    if isa(w,'walk3sagrav')
        animate(w,xs,ts,stepindices,1,stepindices);
    else
        animate(w,xs,stepindices,'frames',anim);
    end
end

if nargout == 0 || doplot
    plotstep(w,xs,ts);
    if nargout == 0, return; end
end	

xd = xdtemp;
