function [h, jntlocs, comlocs] = drawmodel(w,x2,posnoffset2,h,varargin)
% drawmodel(w, x, posnoffset) draws the 2-d walking model,
% with the model offset relative to the origin by posnoffset (x and y)
% h = drawmodel(w) draws the model at its fixed point and a position zero
%   and returns the handle to the graphic
% drawmodel(w, x, posnoffset, h) takes the graphics handle h and redraws the
%   model in that handle (useful for animation)
% drawmodel(w, x, posnoffset, h, 'modellinewidth', 4) applies the named
%   options to the model.
%
% Available options:
%   modellinewidth   thickness of lines for feet and legs, in points
%   footn            number of segments to draw for foot
%   jointsize        size of hip joint circle, in points
%   showgrfs         (0,1) draw arrows for ground reaction force
%   showcomvel       (0,1) draw arrows for COM velocity

% Art Kuo

if nargin < 2 || isempty(x2)
    x2 = get(w, 'xstar'); % assume fixed point if no other info
end

if nargin < 3 || isempty(posnoffset2)
    posnoffset2 = [0; 0]; % if nothing specified, draw it at origin
end

parms = get(w,'parms');
L = parms.L; R = parms.R; C = parms.C; Mp = parms.Mp; M = parms.M;

% drawing parameters
footn = 12; % how many segments to draw in the curved foot
modellinewidth = 4; % how wide are the line segments
jointsize = 6; % how big is the marker at each joint, in points

% arrow parameters
aang = pi/6; scale = 0.05; scale2 = 2; vx2 = 0.4; vy2 = 1.2;

showcomvel = 0; showgrfs = 0;
internal = 0;

% Check for options
if nargin > 4
    property_argin = varargin;
    secondargument = property_argin{1};
    % Step through the optional arguments
    while length(property_argin) >= 2,
        prop = property_argin{1};
        val = property_argin{2};
        property_argin = property_argin(3:end);
        switch prop
            case 'modellinewidth'
                modellinewidth = val;
            case 'footn'
                footn = val;
            case 'jointsize'
                arrowflag = val;
            case 'showgrfs'
                showgrfs = val;
            case 'showcomvel'
                showcomvel = val;
            case 'internal'
                internal = val;
                
        end % switch
    end % looping through arguments
end % if there are optional arguments

for internal_loop = 1:1+internal
    x = x2((1:4)+4*internal_loop-4);
    posnoffset = posnoffset2(:, internal_loop);
    
    % angles of stance leg, swing leg
    q1 = x(1); q2 = x(2);
    
    contactpoint = -q1*R;
    Rot1 = [cos(q1) -sin(q1); sin(q1) cos(q1)];
    Rot2 = [cos(q2) -sin(q2); sin(q2) cos(q2)];
    
    % A foot
    alpha = 0.8; % need to decide on foot length
    % generic foot goes from -alpha to +alpha
    footang = linspace(-alpha, +alpha+0.4, footn);
    foot = R*[sin(footang); -cos(footang)];
    
    % Points of interest, arranged as [x;y] locations
    
    stankle = [contactpoint + posnoffset(1); posnoffset(2) + R]; % ankle is hub of wheel
    stheel = stankle + [R*sin(q1); -R*cos(q1)]; % heel is end of leg at curved foot
    hip = stheel + [-L*sin(q1); L*cos(q1)];   % hip is where legs and trunk meet
    swheel = hip + [L*sin(q2); -L*cos(q2)];   % heel is end of leg at curved foot
    swankle = swheel + [-R*sin(q2); R*cos(q2)]; % ankle is hub of wheel
    stlegc = hip + [(L-C)*sin(q1); -(L-C)*cos(q1)];   % hip is where legs and trunk meet
    swlegc = hip + [(L-C)*sin(q2); -(L-C)*cos(q2)];   % hip is where legs and trunk meet
    
    % Leg segments, arranged as [x row; y row]
    stfoot = Rot1*foot;
    stfoot(1,:) = stfoot(1,:) + stankle(1); stfoot(2,:) = stfoot(2,:) + stankle(2);
    stleg = [stheel hip];
    swleg = [hip swheel];
    swfoot = Rot2*foot;
    swfoot(1,:) = swfoot(1,:) + swankle(1); swfoot(2,:) = swfoot(2,:) + swankle(2);
    
    if showcomvel
        comloc = (Mp*hip + M*stlegc + M*swlegc) / (Mp+2*M);
        comvel = [(-M*(R+(C-R)*cos(q1))-(Mp+M)*(R+(L-R)*cos(q1)))*x(4) + ...
            +M*(C-R)*cos(q2)*x(5) ;
            (-M*(R+(C-R)*sin(q1))-(Mp+M)*(R+(L-R)*sin(q1)))*x(4) + ...
            +M*(C-R)*sin(q2)*x(5)];
    end
    
    if showgrfs
        [grf,grt,cop] = groundreactionforces(w, x, 0);
        grfloc = [contactpoint + posnoffset(1); posnoffset(2)] + [cop;0];
    end
    
    % form matrices of the x and y positions of line segments,
    % first row is x, second is y, and
    % each column is another point on the line
    
    if nargin < 4 || isempty(h) || length(h)<internal_loop*10+1 % no handle exists, so create new lines
        if(internal_loop == 1+internal) % real model
           legColor = [0 0 1];
        else
            legColor = [.6 .6 .8];
        end
        h(internal_loop*10+1) = line(stfoot(1,:),stfoot(2,:),'color',legColor,'LineWidth',modellinewidth,'EraseMode','normal');
        h(internal_loop*10+2) = line(stleg(1,:),stleg(2,:),'color',legColor,'LineWidth',modellinewidth,'EraseMode','background');
        h(internal_loop*10+3) = line(swleg(1,:),swleg(2,:),'color',legColor,'LineWidth',modellinewidth,'EraseMode','background');
        h(internal_loop*10+4) = line(swfoot(1,:),swfoot(2,:),'color',legColor,'LineWidth',modellinewidth,'EraseMode','normal');
        h(internal_loop*10+5) = line(hip(1), hip(2), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', ...
            'none', 'MarkerSize', jointsize, 'Marker', 'o'); % mark the hip
        if showcomvel
            h(internal_loop*10+6) = drawcomvel(comvel, comloc);
        end
        if showgrfs
            h(internal_loop*10+7) = drawgrf(grf, grfloc);
        end        
    else
        set(h(internal_loop*10+1),'Xdata',stfoot(1,:),'Ydata',stfoot(2,:));
        set(h(internal_loop*10+2),'Xdata',stleg(1,:),'Ydata',stleg(2,:));
        set(h(internal_loop*10+3),'Xdata',swleg(1,:),'Ydata',swleg(2,:));
        set(h(internal_loop*10+4),'Xdata',swfoot(1,:),'Ydata',swfoot(2,:));
        set(h(internal_loop*10+5),'Xdata',hip(1),'Ydata',hip(2)); % hip joint
        if showcomvel
            drawcomvel(vel, startpt, h(internal_loop*10+6));
        end
        if showgrfs
            drawgrf(grf, grfloc, h(internal_loop*10+7));
        end
        
    end % internal for
end % if

end % function

function hgrf = drawgrf(grf, endpt, h, varargin)
% Draws a grf arrow, with the vector grf
% and with the arrow drawn ending at position endpt
% Also can update the coordinates in an existing handle for the
% arrow, if supplied. Otherwise it draws it from scratch and
% returns the handle.
%
% Optional arguments:
%   gap (0.02)           gap between start of arrow and its origin
%   headsize (0.05)   size of arrowhead
%   arrowang (pi/8)   angle of arrow relative to body
%   scale (0.5)       determines length of arrow

% optional parameters
gap = 0.02; headsize = 0.05; % headsize determines size of arrowhead
arrowang = pi/8;          % angle of arrow head relative to body
scale = 0.5;     % scale determines length of arrow

if nargin > 3 % we have options
    property_argin = varargin;
    % Step through the optional arguments
    while length(property_argin) >= 2,
        prop = property_argin{1};
        val = property_argin{2};
        property_argin = property_argin(3:end);
        switch prop
            case 'gap'
                numsteps = val;
            case 'scale'
                scale = val;
            case 'arrowang'
                arrowang = val;
            case 'headsize'
                headsize = val;
        end
    end
end

forceang = atan2(grf(2), grf(1));
grf = grf*scale;

gx = -gap*cos(forceang) + endpt(1); gy = -gap*sin(forceang) + endpt(2);

arrowx = gx + [-grf(1) 0 -headsize*cos(forceang+arrowang) NaN 0 -headsize*cos(forceang-arrowang)];
arrowy = gy + [-grf(2) 0 -headsize*sin(forceang+arrowang) NaN 0 -headsize*sin(forceang-arrowang)];

if nargin < 3 || isempty(h),
    hgrf = line(arrowx, arrowy, 'LineWidth',2,'Color','Green');
else
    set(h, 'Xdata', arrowx, 'Ydata', arrowy);
end

end % function hgrf

function hvel = drawcomvel(vel, startpt, h, varargin)
% Draws a com velocity arrow, with the velocity vel
% and with the arrow drawn starting at position startpt
% Also can update the coordinates in an existing handle for the
% arrow, if supplied. Otherwise it draws it from scratch and
% returns the handle.
%
% Optional arguments:
%   gap (0)           gap between start of arrow and its origin (unused)
%   headsize (0.05)   size of arrowhead
%   arrowang (pi/8)   angle of arrow relative to body
%   scale (0.5)       determines length of arrow

% optional parameters
gap = 0; headsize = 0.05; % headsize determines size of arrowhead
arrowang = pi/8;          % angle of arrow head relative to body
scale = 0.5;     % scale determines length of arrow

if nargin > 3 % we have options
    property_argin = varargin;
    % Step through the optional arguments
    while length(property_argin) >= 2,
        prop = property_argin{1};
        val = property_argin{2};
        property_argin = property_argin(3:end);
        switch prop
            case 'gap'
                numsteps = val;
            case 'scale'
                scale = val;
            case 'arrowang'
                arrowang = val;
            case 'headsize'
                headsize = val;
        end
    end
end

velang = atan2(vel(2), vel(1));
vel = vel*scale;

arrowx = startpt(1) + [0 vel(1) vel(1)-headsize*cos(velang+arrowang) NaN vel(1) vel(1)-headsize*cos(velang-arrowang)];
arrowy = startpt(2) + [0 vel(2) vel(2)-headsize*sin(velang+arrowang) NaN vel(2) vel(2)-headsize*sin(velang-arrowang)];

if nargin < 3 || isempty(h),
    hvel = line(arrowx, arrowy, 'LineWidth',2,'Color','Black');
else
    set(h, 'Xdata', arrowx, 'Ydata', arrowy);
end

end % function drawcomvel

function out = draw_arrow(startpoint,endpoint,headsize)
%by Ryan Molecke
% accepts two [x y] coords and one double headsize

v1 = headsize*(startpoint-endpoint)/2.5;

theta = 22.5*pi/180;
theta1 = -1*22.5*pi/180;
rotMatrix = [cos(theta)  -sin(theta) ; sin(theta)  cos(theta)];
rotMatrix1 = [cos(theta1)  -sin(theta1) ; sin(theta1)  cos(theta1)];

v2 = v1*rotMatrix;
v3 = v1*rotMatrix1;
x1 = endpoint;
x2 = x1 + v2;
x3 = x1 + v3;
hold on;
fill([x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],[0 0 0]);     % this fills the arrowhead (black)
plot([startpoint(1) endpoint(1)],[startpoint(2) endpoint(2)],'linewidth',2,'color',[0 0 0]);
end