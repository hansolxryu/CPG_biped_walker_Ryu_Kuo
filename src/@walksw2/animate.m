function animate(w,varargin)
% animate(w) 
%   animates the 2-d dynamic walking model walksw2 described by w.
% animate(w, 'numsteps', N) automatically repeats the simulation N times (default 2).
% animate(w, x, framecounts) animates w with the given states given in x
%   (rows of the state vector), and with framecounts containing a list
%   of the number of frames in each step, e.g. [16 16 16] for 3 steps.
% animate(w, 'arrows', 1) also draws arrows for COM velocity
% animate(w, 'save', 'eps') animates w and saves each frame as an 
%   'eps', 'avi', 'mp4', or 'gif' file (or 1, 2, 3, 4 respectively).
% animate(w, 'frames', Nframes) draws a certain number of frames per step
% Other options: 'treadmill' display moving treadmill or walk overground;
%   'framedelay' (0.05) how long to wait between frames of animation;
%   'groundlinewidth' how thick the ground line should be, in points;
%   Also drawmodel options such as 'modellinewidth', 'showcomvel', and
%   'showgrfs' (show arrows for COM velocity or ground reaction forces);      
%   'floorheight' is the vertical height of floor tiles
%   'jointsize' is the size of circular dots placed at joints (pelvis)
%   'footsize' is the size of circular dots placed at the feet 

% Dependency: This calls drawmodel, a separate model-specific function
%   so that animate can remain relatively unchanged across models 

% allow a keypress or button press in the figure to stop the animation
keypressed = 0;
set(gcf,'KeyPressFcn', @keypress, 'ButtonDownFcn', @keypress); 

treadmill = 1; % by default, draw a treadmill

%parms = get(w,'parms'); % access model parameters
N = get(w, 'N'); L = 1; R = 0;

% dimensions of objects to be drawn
floorheight = 0.02; % default, but can be changed as parameter
modellinewidth = 4; groundlinewidth = 2; jointsize = 10; footsize = 6;

footn = 10; % how many segments to draw in the curved foot

debg = 0; % set to 1 if you want a long pause between frames

numsteps = 2; saveflag = 0; nframes = 16; % default values
framedelay = 0.05;

showcomvel = 0; showgrfs = 0;

x0 = []; framecounts = []; % default is to use fixed point

setFinishLikeStart = NaN;

% Allow for flexible input arguments, particularly allowing second argument
% to be x, a list of state vectors to animate. Otherwise assume the
% inputs are property/value pairs.
if nargin == 0
  error('animate: need a walk object as argument');
elseif nargin > 1 % optional arguments
  % figure out if the second argument is an initial condition vector, or a
  % bunch of rows of states
  property_argin = varargin;
  secondargument = property_argin{1};
  if isa(varargin{1}, 'double')      % second argument appears to be a number
    if length(secondargument) == N   % and it's a vector, meaning an initial condition
      x0 = secondargument;
      property_argin = property_argin(2:end); 
    elseif size(secondargument,1) > 1 && size(secondargument,2) == N % it's a matrix of states
      xs = varargin{2}; framecounts = varargin{3}; numsteps = length(framecounts);
      property_argin = property_argin(3:end);
    else
      error('animate: unknown second argument')
    end
  end
  % Step through the optional arguments
  while length(property_argin) >= 2,
    prop = property_argin{1};
    val = property_argin{2};
    property_argin = property_argin(3:end);
    switch prop
      case 'numsteps'
        numsteps = val;
      case 'save'  % save as 'eps', 'avi', or 'gif'
        if ~isnumeric(val) % a string
          saveflag = val;
        else % a number
          switch val  % legacy allows for number values
            case 1
              saveflag = 'eps';
            case 2
              saveflag = 'avi';
            case 3
              saveflag = 'mp4';
            case 4
              saveflag = 'gif';
          end
        end 
        numsteps = 1; % only save one step by default
      case 'frames'
        nframes = val;
      case 'treadmill'
        treadmill = val;
      case 'showgrfs'
        showgrfs = val;
      case 'showcomvel'
        showcomvel = val;
      case 'framedelay'
        framedelay = val;
      case 'modellinewidth'
        modellinewidth = val;
      case 'groundlinewidth'
        groundlinewidth = val;
      case 'floorheight'
        floorheight = val;
      case 'finishLikeStart'
        setFinishLikeStart = val;
    end
  end
end

% do you want the animation to stop on the first
% frame of a step, so that it looks like it did in the beginning?
% If just viewing on screen, answer is usually yes.
% Or if saving the movie, then usually it's going to be looped, so
% probably end on the frame just before heelstrike, so that there's
% a smooth animation across heelstrike.
if isnan(setFinishLikeStart) % no option set, so let's guess best behavior
  if ischar(saveflag) % saving as any of eps, gif, avi
    finishLikeStart = 0;
  else
    finishLikeStart = 1;
  end
else % a value was given
  finishLikeStart = setFinishLikeStart;
end

% send this to drawmodel
drawmodeloptions = {'modellinewidth', modellinewidth, 'footn', footn, ...
  'jointsize', jointsize, 'footsize', footsize, 'showcomvel', showcomvel, 'showgrfs', showgrfs};

if isempty(framecounts) % we've haven't been supplied with an array of states
  % so do a simulation and figure out how many frames per step
  [xe,te,xs,ts] = onestep(w, x0, 'anim', nframes); framecounts = length(xs);
end

speedinfo = walkspeed(w); sl = speedinfo.steplength;

% Estimate range of walking
distance = numsteps*sl; 

if treadmill % make room for treadmill or overground steps
  xlimit = [-(sl+R)*1.1 sl*1.2];    % horizontal
else
  xlimit = [-(sl+R) distance+0.5*R];% make room for the whole walk
end    
ylimit = [-floorheight-0.04 1.02*L];        % vertical limit
if showgrfs, ylimit(1) = ylimit(1) - 0.1; end; % extra room for grf arrows

% Initialize figure with equal axes with correct limits
clf; set(gcf, 'color', [1 1 1]); 
set(gca,'DataAspectRatio',[1 1 1],'Visible','off','NextPlot','Add','XLim',xlimit,'YLim',ylimit);
set(gca,'Units', 'Points'); 
axesinpoints = get(gca,'Position'); ptstodata = diff(xlimit)/axesinpoints(3);

if saveflag % need to provide a bounding box to make the frames consistent when doing animation
  hboundingbox = rectangle('Position',[xlimit(1) ylimit(1) diff(xlimit) diff(ylimit)],'Visible','On','FaceColor','w','EdgeColor','none');
  switch saveflag % and initiate video file
    case 'avi'
      videofile = VideoWriter(['walkanim.' saveflag]);
      open(videofile)
    case 'mp4'
      videofile = VideoWriter(['walkanim.' saveflag], 'MPEG-4');
      open(videofile)
  end % switch
end % saving

% Drawmodel will plot the foot's ground contact at (0,0) so let's plot
% the ground and tiles underneath that, taking into account the thickness
% of the lines.
%ygndoffset = -0.5*(groundlinewidth+modellinewidth)*ptstodata; % for curved feet
ygndoffset = -0.5*(groundlinewidth+footsize)*ptstodata; % for point feet
ytileoffset = ygndoffset - 0.5*groundlinewidth*ptstodata;
hgnd = line(xlimit,[1 1]*ygndoffset,'color',[0 0 0],'linewidth',groundlinewidth);

% fill the floor with tiles extending horizontally two extra 
% tiles to the right (one black, one white) so that when the treadmill 
% scrolls left, we always show all tiles

tilelen = sl/4; % tile length, equal to an integer fraction of step length
% tilexs is x positions across the screen, while also giving an odd number 
% of vertices horizontally, for an even number of tiles. 

% Make an odd number of vertices to stretch across the floor, for an 
% even number of tiles.
Nc = floor(diff(xlimit)/tilelen) + 1; if mod(Nc,2)==0, Nc = Nc + 1; end;
tilexs = xlimit(1) + (0:Nc-1)*tilelen;
% Tiles are defined by vertex and face matrices
% where the vertices are rows of x-y-z triples, in two groups of rows.
% The first group is the vertices stretching across the screen at ground
% height, and the second group (below the first) is the same set of
% vertices but at floorheight below ground. 
%  1      2      3      4      5
%  *------*      *------*      *    (We end up not drawing the last
%  |-  I -|      |- II -|      |     vertices on the right since they are 
%  *------*      *------*      *     white) 
%  6      7      8      9     10     

tilevm = [tilexs' ones(Nc,1)*ytileoffset zeros(Nc,1);
  tilexs' ones(Nc,1)*ytileoffset-floorheight zeros(Nc,1)];
% define dark faces from the vertices (white faces are not drawn)
tilefm = repmat((1:2:Nc-1)',1,4) + repmat([0 1 Nc+1 Nc],floor(Nc/2),1);
hgndtiles = patch('Vertices', tilevm, 'Faces', tilefm);

% contactpoint with floor can move backwards if we're on
% a treadmill, or it stays in one place. For additional
% steps it has to move forward at each s2s transition.

startpoint = 0; % start us forward so it's approximately
% centered on the treadmill
if treadmill == 0, startpoint = -xs(1,1)+R; end % or try zero

h = drawmodel(w, xs(1,:), [startpoint;0], [], drawmodeloptions{:});

cntr = 1; % absolute frame counter

for j = 1:numsteps+finishLikeStart % Loop through each step
  % finishLikeStart adds +1 step, to end on frame similar to beginning
  if treadmill % we'll move the model with the treadmill
    modelorigin = startpoint; 
    dx = sl / framecounts(min(end,j)); % how much to move the treadmill
  else         % or allow it to advance down the floor
    modelorigin = startpoint + (j-1)*sl;
    dx = 0;
  end
  framesthisstep = framecounts(min(end, j));
  if finishLikeStart && j == numsteps+1
    framesthisstep = 1;
  end
  for i = 1:framesthisstep % each step has this many frames
    drawmodel(w, xs(i,:), [modelorigin; 0], h, drawmodeloptions{:});
    set(hgndtiles,'Vertices', tilevm); % update treadmill vertices
    
    modelorigin = modelorigin - dx; % move model to left    
    tilevm(:,1) = tilevm(:,1) - dx; % and move treadmill to the left
    
    if tilevm(2,1) < xlimit(1), % scrolled a whole tile to left
      tilevm(:,1) = tilevm(:,1) + 2*tilelen;  % reset two tiles to the right
    end
        
    drawnow; pause(framedelay)

    % save animation as eps frames or animated gif
    if saveflag 
      switch saveflag
        case 'eps'
          print('-depsc2',sprintf('walk%02d',cntr));
          % Note use -depsc2 to save in vector format. To load into
          % Adobe Flash, first use Illustrator to import the eps, then
          % export in Flash format.
        case {'avi', 'mp4'}
          frame = getframe;
          writeVideo(videofile, frame);
        case 'gif'
          frame = getframe;   % to save animated gif
          [im,map] = rgb2ind(frame.cdata, 256, 'nodither');
          if cntr == 1 % first frame
            imwrite(im, map, 'walkanim.gif', 'gif', 'Loopcount', Inf);
          else
            imwrite(im, map, 'walkanim.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0);
          end
      end % switch
    end  % if save animation
    
    if debg, pause, end;
      cntr = cntr + 1;
    if keypressed
      return
    end
  
  end % frame loop
  
end % numsteps loop

if strcmp(saveflag, 'avi') || strcmp(saveflag, 'mp4')
  close(videofile); % close video
end

function keypress(src,eventdata)
  keypressed = 1; % in scope of animate function
end
  
end % animate function


