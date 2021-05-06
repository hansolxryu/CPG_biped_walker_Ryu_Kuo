function w = walksw2(varargin)
% WALKSW2 Create a simple 2-d walking gait
%   w = walksw2 creates a default sw2 gait
%   w = walksw2('normal') creates a normal sw2 gait, speed 1.25 m/s, 1.8 Hz
%   w = walksw2('slow') creates a slow sw2 gait, speed 0.1, Kp = 0    
%   w = walksw2('normalgrav') creates a normal gait speed but powered by
%       gravity alone

switch nargin
case 0
  % if no arguments, create a default gait from a good initial guess
  % Here are default parameters for a very slow walk
  w.parms.P = 0.001; w.parms.Kp = 0.125; 
  w.parms.gamma = 0;
  w.N = 4; % 2-d simple walking model has 4 states
  w.xstar = []; % if fixed point not known
  w.xstar = [0.030701 0.061402 -0.0325631 -6.13654e-005];
  w.parent = [];
  w = class(w, 'walksw2');
  [w, cnvrg] = gradsearch(w, [], 'info', 0);
case 1
  if isa(varargin{1}, 'walksw2')
    w = varargin{1}; 
  elseif isempty(varargin{1}) % we want an empty object
    % empty objects are used when another object inherits methods from
    % the walksw2 class
    w.parms = [];
    w.N = []; w.xstar = []; w.parent = [];
    w = class(w, 'walksw2');
  elseif strcmp(varargin{1}, 'normal')
    % Here is a gait that gives equivalent of speed 1.25 m/s, 1.8 Hz
    % P = 0.1877, Kp = 2.5127, xstar = [ 0.355461 0.710921 -0.505623 -0.122482 ]
    % saved in a .mat file as a struct
    fname = which('walksw2/private/normalgait.mat'); loadvar = load(fname);
    w = class(loadvar.wstruct, 'walksw2'); 
  elseif strcmp(varargin{1}, 'normalgrav')
    % Here is a gait that gives equivalent of speed 1.25 m/s, 1.8 Hz
    % with passive gravity power only (no active push-off)
    fname = which('walksw2/private/normalgravgait.mat'); loadvar = load(fname);
    w = class(loadvar.wstruct, 'walksw2'); 
  elseif strcmp(varargin{1}, 'slow')
    % Here is a gait that gives speed 0.1, with
    % no hip spring
    % P = 0.04, Kp = 0, xstar = [ 0.194589 0.389178 -0.202961 -0.0151771 ]
    fname = which('walksw2/private/slowgait.mat'); loadvar = load(fname);
    w = class(loadvar.wstruct, 'walksw2');
  elseif isstruct(varargin{1}) % turn a structure into a class
    w = class(varargin{1}, 'walksw2');
  else
    fprintf(1,[' ''' varargin{1} ''' is not a valid saved object name.  Returning base constructor.']);
    w = walksw2;
  end
end
