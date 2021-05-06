function [bcres, steplength, stepperiod] = bvpbcfun(xa, varargin)
% BVPBCFUN evaluates the residual of the boundary condition for the
% boundary value problem solver, bvpsearch.
% For walksw2, there are five boundary conditions: four from the
% periodicity required after step-to-step transition, and one for
% the double support condition. There is one unknown bvp parameter,
% the step period T. This is separate from the walk object's "known"
% parameters structure parms.
%   bcres = bvpbcfun(xa, xb, T, walk, walkparms) as used by bvp4c
%   bcres = bvpbcfun(xa, xb, T, walk, walkparms) for
%     multipoint problems.
% Alternative calls: bvpbcfun(walk [, sol]) takes a walk object with
% or without solution and returns the residual, as called by the user.
% 
% [bcres, steplength, stepperiod] = bvpbcfun... also returns extra
%   gait info.
% Note that xa and xb are expected to be column vectors. Also note that
% other walk models might have more than one unknown bvp parameter.

% Art Kuo 6/2009

if nargin == 5, % normal call, assign the variables
  [xb, unknownparms, walk, walkparms] = deal(varargin{:});
else
  if isobject(xa) % we're given a walk object as first argument
    walk = xa; % in which case, the second argument may or may not be
    % a solution structure
    if nargin == 1   % we're not given a solution structure
      sol = bvpsetup(walk, get(walk,'xstar')); % so let bvpsetup use its xstar
    elseif nargin == 2 && isstruct(varargin{1}) % we are given solution
      sol = varargin{1};
    else
      error('bvpbcfun called with incorrect input arguments');
    end
    walkparms = get(walk, 'parms');
  end
  % Now we have the solution, so set it up to find boundary condition 
  % residual
  % Look for duplicates in sol.x, indicating region boundaries
  xa = sol.y(:,1);
  xb = sol.y(:,end);
  unknownparms = sol.parameters;
end

T = unknownparms(1); % walksw2, only one unknown parameter for regular
% bvpsearch, but unknownparms may contain other items for bvpgaitspeed

% take the state at the end condition and apply the step-to-step transtion
xnew = s2stransition(xb', walk, walkparms)';
% (Note that s2stransition returns a row vector, so we transpose it)
% we want xnew to match the beginning condition, plus there is one
% additional boundary condition which is the Poincare section, or
% heel strike event
bcres = [xnew - xa; heelstrikeevent(T, xb', walk, walkparms)];

steplength = -2*sin(xb(1,end)); % use ending condition to determine step length
stepperiod = T;

end