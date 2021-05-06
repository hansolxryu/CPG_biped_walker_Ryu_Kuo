function w = walksw2cpg(varargin)
% WALKSW2 Create a simple 2-d walking gait
%   w = walksw2 creates a default sw2 gait
%   w = walksw2('normal') creates a normal sw2 gait, speed 1.25 m/s, 1.8 Hz
%   w = walksw2('slow') creates a slow sw2 gait, speed 0.1, Kp = 0
switch nargin
    case 0
        % if no arguments, create a default gait
        % Here are default parameters
        
        x0 = [0.25 -0.25 -1.0 -0.6];
        wsw2 = set(walksw2, 'xstar', x0,...
            'gamma', 0, 'P', 0, 'Kp', 2,'N',4);
        w.parent = 'walksw2';
        
        w.xstar = x0;
        w.N = 4;
        w.parms.testParm = 0;  %% added for test inheritance
        w.parms.noiseVector = zeros(6,1);
        w.parms.Xr = zeros(4,1);
        w.parms.Lgain = zeros(4,4);        
        
        % anthropomorphic model % 
        w.parms.g = 1;
        w.parms.Mp = 0.68; 
        w.parms.M = 0.16;
        
        w.parms.L = 1; 
        w.parms.C = 0.645; 
        w.parms.R = 0.3;
        w.parms.Ip = 0; 
        w.parms.Il = w.parms.M * (0.326*w.parms.L)^2;
        
        w.parms.Tst = -0.35;
        w.parms.ctrlPortion = 1;
        w.parms.theta0 = 0;
        w.parms.damping = 0;
        w.parms.swingCtrl = 0;
        w.parms.changeST = 0;
        w.parms.testParmMdl = [];
        w.parms.internalState = 1;

        w.parms.A0 = [];
        w.parms.CutY = 0;
        w = class(w, 'walksw2cpg', wsw2); 
%         [w, cnvrg] = gradsearch(w, [], 'info', 0, 'criterion', 1e-9);
        
    case 1
        if isa(varargin{1}, 'walksw2')
            w = varargin{1};
        elseif isempty(varargin{1}) % we want an empty object
            % empty objects are used when another object inherits methods from
            % the walksw2 class
            w.parms = [];
            w.N = []; w.xstar = []; w.parent = [];
            w = class(w, 'walksw2ctrl');
            %   elseif strcmp(varargin{1}, 'normal')
            %     % Here is a gait that gives equivalent of speed 1.25 m/s, 1.8 Hz
            %     % P = 0.1877, Kp = 2.5127, xstar = [ 0.355461 0.710921 -0.505623 -0.122482 ]
            %     fname = which('walksw2/private/normal.mat'); loadvar = load(fname);
            %     w = loadvar.w;
            %   elseif strcmp(varargin{1}, 'normalgrav')
            %     % Here is a gait that gives equivalent of speed 1.25 m/s, 1.8 Hz
            %     % with gravity power only
            %     fname = which('walksw2/private/normalgrav.mat'); loadvar = load(fname);
            %     w = loadvar.w;
            %   elseif strcmp(varargin{1}, 'slow')
            %     % Here is a gait that gives speed 0.1, with
            %     % no hip spring
            %     % P = 0.04, Kp = 0, xstar = [ 0.194589 0.389178 -0.202961 -0.0151771 ]
            %     fname = which('walksw2/slowgait.mat'); loadvar = load(fname);
            %     w = loadvar.w; w = class(w, 'walksw2ctrl');
        else
            fprintf(1,[' ''' varargin{1} ''' is not a valid saved object name.  Returning base constructor.']);
            w = walksw2ctrl;
        end        
end
% cnvrg
end
