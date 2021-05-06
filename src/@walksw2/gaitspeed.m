function [speed, sl, sf, energy] = gaitspeed(w, varargin)
% [speed, steplength,stepfreq] = gaitspeed(w) returns the speed, steplength,
%   and step frequency in dimensionless units for a gait w.
%   w is a simplest walking model in 2-d, walksw2
%
% [speed, steplength, stepfreq, energy] = gaitspeed(w) also adds some info on
%   energy expenditure, in the form of a structure with
%   energy.cot = mechanical cost of transport (energy per weight and distance)
%   energy.power = average mechanical power
%     Where this is computed only if the energy output is requested.
%
% [speed, steplength, stepfreq] = gaitspeed(w, xstar) uses the initial
%   condition xstar, instead of looking up the one stored in w.
%
% [speed, steplength, stepfreq] = gaitspeed(w, xstar, xend, tend)
%   optinally provide the end state and time to eliminate the need to do
%   a simulation.

% Modifications: added energetics as an extra output (Art, 7/24/2008)
%   But this should be considered necessary for all models
% 
%   Added a shortcut to allow xstar, xend, tend to be provided to
%     reduce need for simulations (Art, 10/2010)

xs = get(w, 'xstar'); gamma = get(w, 'gamma');

if nargin > 1 && ~isempty(varargin{1}) % use the provided initial condition if it exists
  xs = varargin{1};
end
  
if nargin > 1 && length(varargin) == 3 % end state and time provided
  xe = varargin{2}; te = varargin{3};
else % need to do a simulation
  [xe,te,x,t,e] = onestep(w);
end

sl = 2*sin(xs(1)); % these apply to the simplest walker only
sp = sl / te;
sf = 1/te;

if nargout == 0
    fprintf(1, ['speed: ' num2str(sp) '   step length: ' num2str(sl)  '   step frequency: ' num2str(sf) '\n'])
    return
end

speed = sp;

if nargout == 4 % only compute energy if requested as output
  energy.cot = e.pushoffwork / sl + sin(gamma);
  energy.power = energy.cot * speed; % average rate of COM work
end
