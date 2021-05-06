% Contents of walking class library
%
% TRUNK/
% demowalking           runs several demonstrations of walksw2 and walk2 gaits
% extrapolategait       library function for guess new fixed points given
%                       existing ones, useful for parameter studies
% interpolategait       library function for producing intermediate gaits
%                       between two extremes, useful for finding fixed points
% findgait              given a known gait and fixed point, try to find another
%                       gait with different parameters and unknown fixed point
% findgaitspeed         find a gait that matches speed and step length
%                       requirements
% findgaitspeeddutyfactor  find a gait that matches speed, step length, 
%                       and duty factor requirements (for springy models)
% gradient              return partial derivative of a gait about fixed point
% gradsearch            perform first-order search for fixed point (usually
%                       preferrable to use findgait)
% manystep              runs onestep function over many steps and returns
%                       results graphical and numerical results
% parmstudy1d           performs one-dimensional parameter study on a range of
%                       specified parameter values
% stability             return eigenvalues and eigenvectors from fixed pt gradient
% stepsearch            given a set of initial conditions, try to find a fixed
%                       point by taking many steps until convergence
% walkspeed             returns the speed, step length, step frequency, and
%                       optionally duty factor of a gait
%
% @walksw2/        simplest planar walking model       
%   @walk2/        anthropomorphic 2-d model
%     @walk2k/     anthropomorphic 2-d model with knees
%     @walk2off/   anthropomorphic 2-d model with forward offset of the foot
%     @walk3/      anthropomorphic 3-d model (not yet completed)
%   @walk2e/       anthropomorphic 2-d model with numerical derivation of EOM
% @hybridsys/      constructor class for hybrid simulations
%   @walk2h/         hybrid simulation of anthropomorphic 2-d model
%   @walk2kh/        hybrid simulation of anthropomorphic 2-d model with knees

% DOCUMENTATION/        contains Word documentation files