function [varargout] = walkspeed(w, varargin)
% infostructure = walkspeed(w) returns the speed, step length,
%   step frequency, and optionally duty facor in dimensionless 
%   units for a gait w. Output is in the form of a structure with 
%   fields 'speed', 'steplength', 'stepfreq', and optional 'dutyfact'
% walkspeed(w) with no output arguments prints out the speed, step length,
%   and step frequency
% walkspeed(w, 'plotfigure', 1) plots the angles for a step and then prints the
%   speed, step length, and step frequency on the figure
% [speed, steplength, stepfreq] = walkspeed(w) or
% [speed, steplength, stepfreq, dutyfact] = walkspeed(w) returns the same 
%   information as scalars, if the user does not want a structure

% Art Kuo  v 1.0  (04/07/2004)

plotfigure = 0;

opt_argin = varargin;         % parse the options
onestepoptions = {};
while length(opt_argin) >= 2,
 opt = opt_argin{1};
 val = opt_argin{2};
 opt_argin = opt_argin(3:end);
 switch opt
   case 'plotfigure'
     plotfigure = val;
   case {'AbsTol','RelTol', 'tf'} % options for onestep
     onestepoptions = {onestepoptions{:}, opt, val};
   otherwise
     warning(['walkspeed options: plotfigure']);
 end
end

try
    [speed,steplength,stepfreq,dutyfact] = gaitspeed(w);
catch
   errmsg = lasterr;
   if(strfind(errmsg, 'Too many output arguments.'))
      [speed,steplength,stepfreq] = gaitspeed(w);
      dutyfact = 0.5;
   else
      rethrow(lasterror);
   end
end

if plotfigure % do a onestep, plot the angles, and show the speed info
 onestep(w, [], onestepoptions);
 if dutyfact == 0.5
 outstring = sprintf('speed     = %4.2f\nstep len  = %4.2f\nstep freq = %4.2f', speed, steplength, stepfreq);
 else
 outstring = sprintf('speed     = %4.2f\nstep len  = %4.2f\nstep freq = %4.2f\nduty fact = %4.2f', speed, steplength, stepfreq, dutyfact);   
 end
 xl = get(gca, 'XLim'); yl = get(gca, 'YLim');
 % try to put the info on the upper left corner
 x = xl(1) + 0.15*(xl(2)-xl(1)); y = yl(2) - 0.03*(yl(2)-yl(1));
 text(x, y, outstring, 'VerticalAlignment', 'top');
end

sinfo.speed = speed;
sinfo.steplength = steplength;
sinfo.stepfreq = stepfreq;
sinfo.dutyfact = dutyfact;
varargout{1} = sinfo;

if nargout >= 3 % user apparently wants three scalar outputs instead of a structure
 varargout{1} = speed;
 varargout{2} = steplength;
 varargout{3} = stepfreq;
 varargout{4} = dutyfact;
end