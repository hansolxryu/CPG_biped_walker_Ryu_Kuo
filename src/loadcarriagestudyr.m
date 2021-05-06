%% Load carriage study with arc feet
% Using walksw2dhlcr, the double-humped simplest walker
% with load carriage.
%
% Walksw2dhlcr contains a parameter M for body mass and an arc foot R = 0.3,
% and computes all forces
% and energies including this parameter. It is thus necessary to scale the
% total energy E according to total mass. We also perform a set of 
% simulations where speed and step length are kept constant, and E and Kp
% are scaled. These are referred to as ScaledE and ConsV, respectively.
%
% The main finding is that both models yield approximately constant
% speed, but in different ways. ScaleE tends to take longer steps that
% compensates for longer double support/slower step frequencies. ConsV
% needs increased Kp to go faster, but also tends to increase E. Since
% both need about the same energy, it seems that ScaledE is easier.
%
% The effect of added mass is to cause double support time to increase.
% COM work increases for most phases as well, slightly faster than
% linearly. 

% Start with a normal walksw2dh
% let's reduce the step length slightly, to match our data a bit better
wdhlcr = walksw2dhlcr(walksw2dh('normal'));
%wnew3 = findgait(wdhlcr, set(wdhlcr,'R',0.3)); R = 0.3;  % the slow way
wnew3 = gradsearch(set(wdhlcr,'R',0.3,'Kp',7.12077,'E',1.04419','Kleg',18.3053, ...
  'xstar',[ 0.443585 0 -0.229659 0.0506988 -0.3547 0.195584 -0.445795 0.00431439 ]));
wnew3 = findgaitspeed(wnew3, 'speed', 0.4, 'steplength', 0.65, 'parmvary1', 'E', 'parmvary2', 'Kp')
wstar3 = wnew3; R = 0.3;

walksw2dhlcrDirectory = fileparts(which('walksw2dhlcr'));

%% Find some higher stiffness gaits
% Here's one with R = 0.3, Kleg = 35
%try1=findgait(wstar3,set(wstar3,'Kleg',35)); % the slow way
try1 = gradsearch(set(wstar3,'xstar',[ 0.397562 0 -0.269338 0.0304558 -0.394432 0.158509 -0.447525 -0.0419979 ],...
  'Kleg',35,'Kp',4.65672,'E',1.0567,'steplength',0.649251,'R',0.3));
wstar35=findgaitspeed(try1,'speed',0.4, 'steplength',0.65,'parmvary1','E','parmvary2','Kp')
% Or you can also use the boundary value solver:
wstar35b=findgaitspeedb(try1,'speed',0.4, 'steplength',0.65,'parmvary1','E','parmvary2','Kp')

%% Parameter study increasing mass, scaling E with mass
% This turns out to cause step length to increase slightly, but
% overall speed stays about the same because step frequency goes down a
% bit. This is not too different from what happens in experiment.
Ms = 1:0.05:1.4; w = wstar35; 
Estar = get(w,'E');
if exist([walksw2dhlcrDirectory '/masses.mat']) % pre-saved results
  load([walksw2dhlcrDirectory '/masses.mat'])
else % do the study
clear masses
for i=1:length(Ms)
  [w,cnvrg] = findgait(w, set(w, 'Mp', Ms(i), 'E', Estar*Ms(i)));
  masses.w(i) = w;
  masses.cnvrg(i) = cnvrg;
  if cnvrg
    speedinfo = walkspeed(w);
    masses.speed(i) = speedinfo.speed;
    masses.steplength(i) = speedinfo.steplength;
    masses.stepfreq(i) = speedinfo.stepfreq;
    masses.Ms(i) = Ms(i);
    masses.energies(i) = energies(w);
  end
end

% look at energy over a single step to verify energy conservation
clf; hold on;
for i = 1:length(Ms)
  energies(masses.w(i));
end

% post-processing to determine double support time, com work over phases,
% and other info
clf; hold on;
for i = 1:length(Ms)
  masses.Kp(i) = get(masses.w(i),'Kp');
  masses.E(i) = get(masses.w(i), 'E');
  [xc,tc,xs,ts,es,ies] = onestep(masses.w(i));
  toTime = find(ies == 1, 1);     % toe-off time (event 1 for walksw2dhlcr)
  masses.ds(i) = ts(toTime) / tc; % double support as fraction of stride
  masses.dsAbs(i) = ts(toTime);   % double support as an absolute
  [grf,grt,cop,comworkrate,finfo] = groundreactionforces(ts, xs, masses.w(i));
  masses.comphases(i,:) = finfo.comphases; % CO, RB, PL, PO in array
  masses.totalwork(i) = finfo.totalwork;   % total negative work
  masses.dswork(i,:) = finfo.dswork;         % double support work 
  masses.minmaxfy(i,:) = finfo.minmax;    % includes min and max vert forces
  ycom = R+cos(xs(:,1)).*(1-R-xs(:,2));
  masses.comydisp(i) = max(ycom)-min(ycom); % vert COM displacement
  masses.delta(i) = finfo.delta;
  masses.vminus(i) = finfo.vminus;
  plot(ts, comworkrate);
end
% save([walksw2dhlcrDirectory '/masses.mat'], 'masses')
clf;plot(masses.Ms, masses.comphases) % show that work per phase increases with mass
xlabel('Total mass'); ylabel('COM work'); legend('CO','RB','PL','PO');
end % do the study

%% Parameter study making speed and step length match empirical
% Now do a study varying mass while fitting step length and frequency to
% the empirical data.
% This is found by varying total energy E, and hip spring Kp to
% maintain the empirical gait. This requires a slightly increasing Kp 
if exist([walksw2dhlcrDirectory '/massesemp.mat']) % pre-saved results
  load([walksw2dhlcrDirectory '/massesemp.mat']);
else % do the study
speed = 0.4; % keep the same standard speed
% here's the slope and offset of step length, where s = 0.65 at M = 1, and
% s = 0.7 at M = 1.5, taken from our load carriage data.
lengtheffect = [(0.7-0.65)/(1.5-1) 0.65-(0.7-0.65)/(1.5-1)];
lengths = Ms*lengtheffect(1) + lengtheffect(2);
clear massesemp % varying masses, constant step length and frequency
for i=1:length(Ms)
  [w,cnvrg] = findgaitspeed(masses.w(i), speed, lengths(i), 'parmvary1', 'E', 'parmvary2', 'Kp');
  massesemp.w(i) = w;
  massesemp.cnvrg(i) = cnvrg;
  if cnvrg
    speedinfo = walkspeed(w);
    massesemp.speed(i) = speedinfo.speed;
    massesemp.steplength(i) = speedinfo.steplength;
    massesemp.stepfreq(i) = speedinfo.stepfreq;
    massesemp.Ms(i) = Ms(i);
    massesemp.energies(i) = energies(w);
  end
end
clf; hold on;
for i = 1:length(Ms)
  massesemp.Kp(i) = get(massesemp.w(i),'Kp');
  massesemp.E(i) = get(massesemp.w(i), 'E');
  [xc,tc,xs,ts,es,ies] = onestep(massesemp.w(i));
  toTime = find(ies == 1, 1);        % toe-off time
  massesemp.ds(i) = ts(toTime) / tc; % double support as fraction of stride
  massesemp.dsAbs(i) = ts(toTime);   % double support as an absolute
  [grf,grt,cop,comworkrate,finfo] = groundreactionforces(ts, xs, massesemp.w(i));
  massesemp.comphases(i,:) = finfo.comphases; % CO, RB, PL, PO in array
  massesemp.totalwork(i) = finfo.totalwork;   % total negative work
  massesemp.dswork(i,:) = finfo.dswork;         % double support work 
  massesemp.minmaxfy(i,:) = finfo.minmax;    % includes min and max vert forces
  ycom = R+cos(xs(:,1)).*(1-R-xs(:,2));
  massesemp.comydisp(i) = max(ycom)-min(ycom); % vert COM displacement
  massesemp.delta(i) = finfo.delta;
  massesemp.vminus(i) = finfo.vminus;
  plot(ts, comworkrate);
end
clf;plot(massesemp.Ms, massesemp.comphases) % show work increasing per phase
clf;plot(massesemp.Ms, [massesemp.ds; massesemp.dsAbs]) % double support increases with mass
% save([walksw2dhlcrDirectory '/massesemp'], 'massesemp');
end % do the study
%% ScaleE Results: COM motion, ground reaction forces, COM work rate, and
% minimum and maximum vertical forces.
% Increasing mass causes more COM motion, more excursion of vertical
% forces, more peaks in COM work rate, and a taller hodograph
% without much change in horizontal COM velocity.

% plot vs time
clf; 
for i = 1:length(Ms)
  subplot(221); hold on; xlabel('COM x');  ylabel('COM y'); title('COM motion (scaleE)');
  [xc,tc,xs,ts] = onestep(masses.w(i)); 
  px = -sin(xs(:,1)).*(1-R-xs(:,2))-R*xs(:,1);
  py = cos(xs(:,1)).*(1-R-xs(:,2))+R;
  plot(px(1:end-1),py(1:end-1));
  subplot(222); hold on; xlabel('time'); ylabel('Forces'); title('Forces (scaleE)');
  [grf,grt,cop,comworkrate,finfo,hodograph] = groundreactionforces(ts, xs, masses.w(i));
  plot(ts, grf);   
  subplot(223); hold on; xlabel('Total mass'); ylabel('COM work rate'); title('COM work rate (scaleE)');
  plot(ts, comworkrate); 
  subplot(224); hold on; plot(hodograph(:,1), hodograph(:,2)); xlabel('COM x vel'); ylabel('COm y vel'); title('COM Hodograph (ScaleE)');
  axis equal
end
% plot vs gait cycle
% clf
% for i = 1:length(Ms)
%   subplot(221); hold on; xlabel('COM x');  ylabel('COM y'); title('COM motion (ScaleE)');
%   [xc,tc,xs,ts] = onestep(masses.w(i)); plot(-sin(xs(1:end-1,1)).*(1-R-xs(1:end-1,2))-R*xs(1:end-1,1), cos(xs(1:end-1,1)).*(1-R-xs(1:end-1,2))+R);
%   subplot(222); hold on; xlabel('time (stride)'); ylabel('Forces'); title('Forces (ScaleE)');
%   [grf,grt,cop,comworkrate,finfo,hodograph] = groundreactionforces(ts, xs, masses.w(i));
%   plot(0.5*ts/ts(end), grf);   
%   subplot(223); hold on; xlabel('time (stride)'); ylabel('COM work rate'); title('COM work rate (ScaleE)');
%   plot(0.5*ts/ts(end), comworkrate); 
%   subplot(224); hold on; plot(hodograph(:,1), hodograph(:,2)); xlabel('COM x vel'); ylabel('COm y vel'); title('COM Hodograph (ScaleE)');
%   axis equal
% end

%% Empirical-based Results: COM motion, ground reaction forces, COM work rate, and
% minimum and maximum vertical forces.
% Increasing mass causes more COM motion, more excursion of vertical
% forces, more peaks in COM work rate, and a taller hodograph
% without much change in horizontal COM velocity.

% plot vs time
clf; 
for i = 1:length(Ms)
  subplot(221); hold on; xlabel('COM x');  ylabel('COM y'); title('COM motion (EmpSL)');
  [xc,tc,xs,ts] = onestep(massesemp.w(i)); 
  px = -sin(xs(:,1)).*(1-R-xs(:,2))-R*xs(:,1);
  py = cos(xs(:,1)).*(1-R-xs(:,2))+R;
  plot(px(1:end-1),py(1:end-1));
  subplot(222); hold on; xlabel('time'); ylabel('Forces'); title('Forces (EmpSL)');
  [grf,grt,cop,comworkrate,finfo,hodograph] = groundreactionforces(ts, xs, massesemp.w(i));
  plot(ts, grf);   
  subplot(223); hold on; xlabel('Total mass'); ylabel('COM work rate'); title('COM work rate (EmpSL)');
  plot(ts, comworkrate); 
  subplot(224); hold on; plot(hodograph(:,1), hodograph(:,2)); xlabel('COM x vel'); ylabel('COm y vel'); title('COM Hodograph (EmpSL)');
  axis equal
end
% plot vs gait cycle
% clf
% for i = 1:length(Ms)
%   subplot(221); hold on; xlabel('COM x');  ylabel('COM y'); title('COM motion (EmpSL)');
%   [xc,tc,xs,ts] = onestep(massesemp.w(i)); plot(-sin(xs(1:end-1,1)).*(1-R-xs(1:end-1,2))-R*xs(1:end-1,1), cos(xs(1:end-1,1)).*(1-R-xs(1:end-1,2))+R);
%   subplot(222); hold on; xlabel('time (stride)'); ylabel('Forces'); title('Forces (EmpSL)');
%   [grf,grt,cop,comworkrate,finfo,hodograph] = groundreactionforces(ts, xs, massesemp.w(i));
%   plot(0.5*ts/ts(end), grf);   
%   subplot(223); hold on; xlabel('time (stride)'); ylabel('COM work rate'); title('COM work rate (EmpSL)');
%   plot(0.5*ts/ts(end), comworkrate); 
%   subplot(224); hold on; plot(hodograph(:,1), hodograph(:,2)); xlabel('COM x vel'); ylabel('COm y vel'); title('COM Hodograph (EmpSL)');
%   axis equal
% end
%% Empirical-based Results: Speed & step length, energy, Kp, COM work per phase,
% double support time. 
% We find that the two different ways of walking have relatively subtle
% differences. When scaling E with mass, the step length increases somewhat
% and step frequency decreases, although total speed is fairly constant.
% The total E increases roughly with mass. But when keeping speed and
% step length constant, it is necessary to increase Kp. This doesn't have
% much effect on E, so by our measures we would consider it more costly
% not to increase step length. In terms of COM work, there are very similar
% trends except for PO not increasing that much with constant speed. In
% both cases, double support time increases, more with the constant speed
% model. 
clf
subplot(221); plot(masses.Ms, [masses.speed; masses.steplength], '-'); hold on
%plot(massescon.Ms, [massescon.speed; massescon.steplength], '--');
plot(massesemp.Ms, [massesemp.speed; massesemp.steplength], '--');
legend('ScaleE speed', 'ScaleE steplength', 'EmpSL speed', 'EmpSL steplength','Location','South');
xlabel('Total mass'); ylabel('Speed steplength'); title('Speed & steplength');
subplot(222); plot(masses.Ms, [masses.E; masses.Kp]); hold on
%plot(massescon.Ms, [massescon.E; massescon.Kp],'--');
plot(massesemp.Ms, [massesemp.E; massesemp.Kp],'--');
%legend('ScaleE E', 'ScaleE Kp', 'ConsV E', 'ConsV Kp','Location','South');
legend('ScaleE E', 'ScaleE Kp', 'EmpSL E', 'EmpSL Kp','Location','South');
xlabel('Total mass'); title('Energy and Kp');
subplot(223); plot(masses.Ms, [masses.totalwork' masses.comphases]); hold on;
%plot(massescon.Ms, [massescon.totalwork' massescon.comphases], '--'); 
plot(massesemp.Ms, [massesemp.totalwork' massesemp.comphases], '--'); 
legend('total', 'CO', 'RB', 'PL', 'PO','Location','Best');
xlabel('Total mass'); ylabel('Work/step'); title('COM work per phase');
subplot(224); plot(masses.Ms, [masses.ds; masses.dsAbs], '-'); hold on;
%plot(massescon.Ms, [massescon.ds; massescon.dsAbs], '--'); 
plot(massesemp.Ms, [massesemp.ds; massesemp.dsAbs], '--'); 
%legend('ScaleE DS frac', 'ScaleE DS abs', 'ConsV DS frac', 'ConsV DS abs','Location','Best');
legend('ScaleE DS frac', 'ScaleE DS abs', 'EmpSL DS frac', 'EmpSL DS abs','Location','Best');
xlabel('Total mass'); ylabel('Double support'); title('Double support time');
% it seems that humans do swing the leg harder just as they increase their
% overall energy. 


%% Empirical Results: Min & max vertical forces, v^2 delta^2 predictor
% We find that vertical force excursion increases somewhat linearly
% with mass, and the M v^2 delta^2 predictor is pretty decent with mass.
% Note that comapred to walksw2dhlc, I took out the 1/2 that multiplied
% v^2 delta^2.
clf;
subplot(221)
plot(masses.Ms, [masses.minmaxfy diff(masses.minmaxfy,1,2)])
hold on; plot(massesemp.Ms, [massesemp.minmaxfy diff(massesemp.minmaxfy,1,2)], '--');
xlabel('Total mass'); ylabel('Force'); title('Min & max vert forces');
legend('min','max','diff','Location','East');
subplot(222); 
%plot(masses.Ms, masses.comydisp); hold on; % plot COM displacement
%plot(massescon.Ms, massescon.comydisp, '--'); set(gca,'ylim',[0 Inf]);
%xlabel('Total mass'); ylabel('Displacement'); title('Vert COM displ');
plot(masses.Ms .* masses.vminus.^2 .* masses.delta.^2, masses.totalwork); hold on
plot(massesemp.Ms .* massesemp.vminus.^2 .* massesemp.delta.^2, massesemp.totalwork, '--')
plot([0 0.2],[0 0.2],':'); axis equal; set(gca,'Xlim',[0 Inf],'YLim',[0 Inf]);
legend('ScaleE', 'EmpSL','Location','Best');
xlabel('M v^2 delta^2'); ylabel('Total work');
subplot(223);
plot(masses.Ms, masses.delta*180/pi); set(gca,'ylim',[0 Inf]); hold on;
plot(massesemp.Ms, massesemp.delta*180/pi, '--');
xlabel('Total mass'); ylabel('Delta (deg)'); title('Delta vcom');
subplot(224);
plot(masses.Ms, masses.vminus); set(gca,'ylim',[0 Inf]); hold on;
plot(massesemp.Ms, massesemp.vminus, '--'); 
xlabel('Total mass'); ylabel('vminus'); title('vminus & predictor');
plot(masses.Ms, masses.Ms .* masses.vminus.^2 .* masses.delta.^2,'r'); hold on
plot(masses.Ms, massesemp.Ms .* massesemp.vminus.^2 .* massesemp.delta.^2,'r--')
legend('ScaleE v-','EmpSL v-', 'ScaleE Mv^2\delta^2','EmpSL Mv^2\delta^2','Location','West')


%% Study effect of Kleg, keeping the same speed and step length
if exist([walksw2dhlcrDirectory '/ks.mat']) % pre-saved results
  load([walksw2dhlcrDirectory '/ks.mat']);
else % do the study
Kleg = get(wstar3,'Kleg');
% find a smaller Kleg with same gait speed
try1=findgait(wstar3,set(wstar3,'Kleg',15));
lowerstiff=findgaitspeed(try1,'speed',0.4, 'steplength',0.65,'parmvary1','E','parmvary2','Kp')
% also tried some higher stiffnesses and found Kleg = 45 to be pretty 
% close to skipping
%try8 = findgait(wstar3, set(wstar3,'Kleg',44));
%higherstiff = findgaitspeed(try8,'speed',0.4,'steplength',0.65,'parmvary1','E','parmvary2','Kp');
%try9 = findgait(higherstiff, set(higherstiff,'Kleg',45));
%higherstiff2 = findgaitspeed(try9,'speed',0.4,'steplength',0.65,'parmvary1','E','parmvary2','Kp');

Klegs = linspace(15, 43, 15);
w = lowerstiff; clear ks; clf; hold on;
for i = 1:length(Klegs)
  fprintf(1,'i = %d ', i);
  [trialw,cnvrg] = findgait(w, set(w,'Kleg', Klegs(i)));
  if cnvrg
    [w, cnvrg] = findgaitspeed(trialw,'speed',0.4,'steplength',0.65, 'parmvary1','E','parmvary2','Kp');
    if cnvrg
      ks.w(i) = w;
      ks.cnvrg(i) = cnvrg;
      onestep(w);
    else
      break
    end
  else
    break
  end
end

% look at energy over a single step to verify energy conservation
clf; hold on;
for i = 1:length(Klegs)
  energies(ks.w(i));
end

% post-processing to determine double support time, com work over phases,
% and other info
clf; hold on;
for i = 1:length(Klegs)
  ks.Kleg(i) = Klegs(i);
  speedinfo = walkspeed(ks.w(i));
  ks.speed(i) = speedinfo.speed;
  ks.steplength(i) = speedinfo.steplength;
  ks.stepfreq(i) = speedinfo.stepfreq;
  ks.energies(i) = energies(ks.w(i));
  ks.Kp(i) = get(ks.w(i),'Kp');
  ks.E(i) = get(ks.w(i), 'E');
  [xc,tc,xs,ts,es,ies] = onestep(ks.w(i));
  toTime = find(ies == 1, 1);     % toe-off time (event 1 for walksw2dhlcr)
  ks.ds(i) = ts(toTime) / tc; % double support as fraction of stride
  ks.dsAbs(i) = ts(toTime);   % double support as an absolute
  [grf,grt,cop,comworkrate,finfo,hodograph] = groundreactionforces(ts, xs, ks.w(i));
  ks.comphases(i,:) = finfo.comphases; % CO, RB, PL, PO in array
  ks.totalwork(i) = finfo.totalwork;   % total negative work
  ks.dswork(i,:) = finfo.dswork;         % double support work 
  ks.minmaxfy(i,:) = finfo.minmax;    % includes min and max vert forces
  ycom = R+cos(xs(:,1)).*(1-R-xs(:,2));
  ks.comydisp(i) = max(ycom)-min(ycom); % vert COM displacement
  ks.delta(i) = finfo.delta;
  ks.vminus(i) = finfo.vminus;
  %plot(ts, comworkrate);
  plot(ts, grf);
end % for loop
end % do the study
% save([walksw2dhlcrDirectory '/ks.mat'], 'ks')

% Plot trajectories, grfs, com work rate, hodographs
% trajectories
clf; subplot(221); hold all; xlabel('time'); ylabel('qs'); title('Trajectories');
for i = 1:length(ks.Kleg)
  onestep(ks.w(i));
end
% grfs
subplot(222); hold all; xlabel('time'); ylabel('forces'); title('GRFs');
for i = 1:length(ks.Kleg)
  groundreactionforces(ks.w(i));
end
% com work rate
subplot(223); hold all; xlabel('time'); ylabel('work rate'); title('COM work rate');
for i = 1:length(ks.Kleg)
  [grf,grt,cop,comworkrate,finfo,hodograph,t]=groundreactionforces(ks.w(i));
  plot(t,comworkrate);
end
% hodograph
subplot(224); hold all; xlabel('vx'); ylabel('vy'); title('COM hodographs');
for i = 1:length(ks.Kleg)
  [grf,grt,cop,comworkrate,finfo,hodograph,t]=groundreactionforces(ks.w(i));
  plot(hodograph(:,1),hodograph(:,2));
end
axis equal

%% Kleg study: Plot com phases, work, etc. for different Klegs
clf; subplot(221); hold on;
% higher stiffness leads to more rebound and preload
h1=plot(ks.Kleg, ks.comphases); % lower stiffnesses cause rebound and preload to disappear
h2=plot(ks.Kleg, ks.dswork,':'); % higher stiffness means less dswork
h3=plot(ks.Kleg, ks.totalwork,'k--'); % but it also means more total work
xlabel('Kleg'); ylabel('Work'); title('Work vs Kleg');
legend([h1(1); h2(1); h3],'CRPP', 'ds+', 'total');
subplot(222);
plot(ks.Kleg, ks.ds)        % shorter double support
set(gca,'Ylim',[0 Inf]);
xlabel('Kleg'); ylabel('Time (fraction)'); title('Double support time');
subplot(223);
plot(ks.Kleg, ks.Kp, ks.Kleg, 10*ks.E) % Kp goes down, E stays about the same
set(gca,'Ylim',[0 Inf]);
xlabel('Kleg'); ylabel('K and E'); title('Kp and E vs Kleg');
legend('Kp', '10*E')
subplot(224);
plotyy(ks.Kleg, ks.delta, ks.Kleg, ks.vminus);
xlabel('Kleg'); title('deltav and vminus');
legend('deltav','vminus')
% vminus doesn't change much, but deltav increases a lot with stiffness

%% ConsV: Parameter study keeping speed and step length conserved
% Now do a study varying mass while keeping step length and frequency the
% same. This is found by varying total energy E, and hip spring Kp to
% maintain the nominal gait. This requires a slightly increasing Kp.
% In the end, this doesn't seem to fit quite as well as just scaling E
% Deprecated: If you want to run it, set the if statement to 1:
if 0, % set to 1 if you want to see ConsV results
if exist([walksw2dhlcrDirectory '/massescon.mat']) % pre-saved results
  load([walksw2dhlcrDirectory '/massescon.mat']);
else % do the study
speedinfo = walkspeed(masses.w(1)); % keep the same speed as the first gait of masses
clear massescon % varying masses, constant step length and frequency
for i=1:length(Ms)
  [w,cnvrg] = findgaitspeed(masses.w(i), speedinfo.speed, speedinfo.steplength, 'parmvary1', 'E', 'parmvary2', 'Kp');
  massescon.w(i) = w;
  massescon.cnvrg(i) = cnvrg;
  if cnvrg
    speedinfo = walkspeed(w);
    massescon.speed(i) = speedinfo.speed;
    massescon.steplength(i) = speedinfo.steplength;
    massescon.stepfreq(i) = speedinfo.stepfreq;
    massescon.Ms(i) = Ms(i);
    massescon.energies(i) = energies(w);
  end
end
clf; hold on;
for i = 1:length(Ms)
  massescon.Kp(i) = get(massescon.w(i),'Kp');
  massescon.E(i) = get(massescon.w(i), 'E');
  [xc,tc,xs,ts,es,ies] = onestep(massescon.w(i));
  toTime = find(ies == 1, 1);        % toe-off time (event 1 in walksw2dhlcr)
  massescon.ds(i) = ts(toTime) / tc; % double support as fraction of stride
  massescon.dsAbs(i) = ts(toTime);   % double support as an absolute
  [grf,grt,cop,comworkrate,finfo] = groundreactionforces(ts, xs, masses.w(i));
  massescon.comphases(i,:) = finfo.comphases; % CO, RB, PL, PO in array
  massescon.totalwork(i) = finfo.totalwork;   % total negative work
  massescon.dswork(i,:) = finfo.dswork;         % double support work 
  massescon.minmaxfy(i,:) = finfo.minmax';    % includes min and max vert forces
  ycom = R+cos(xs(:,1)).*(1-R-xs(:,2));
  massescon.comydisp(i) = max(ycom)-min(ycom); % vert COM displacement
  massescon.delta(i) = finfo.delta;
  massescon.vminus(i) = finfo.vminus;
  plot(ts, comworkrate);
end
%plot(massescon.Ms, massescon.comphases) % show work increasing per phase
%plot(massescon.Ms, [massescon.ds; massescon.dsAbs]) % double support increases with mass
% save([walksw2dhlcrDirectory '/massescon.mat'], 'massescon')
end % do the study
%% ConsV Results: Speed & step length, energy, Kp, COM work per phase,
% double support time. 
% We find that the two different ways of walking have relatively subtle
% differences. When scaling E with mass, the step length increases somewhat
% and step frequency decreases, although total speed is fairly constant.
% The total E increases roughly with mass. But when keeping speed and
% step length constant, it is necessary to increase Kp. This doesn't have
% much effect on E, so by our measures we would consider it more costly
% not to increase step length. In terms of COM work, there are very similar
% trends except for PO not increasing that much with constant speed. In
% both cases, double support time increases, more with the constant speed
% model. 
clf
subplot(221); plot(masses.Ms, [masses.speed; masses.steplength], '-'); hold on
plot(massescon.Ms, [massescon.speed; massescon.steplength], '--');
legend('ScaleE speed', 'ScaleE steplength', 'ConsV speed', 'ConsV steplength','Location','South');
xlabel('Total mass'); ylabel('Speed steplength'); title('Speed & steplength');
subplot(222); plot(masses.Ms, [masses.E; masses.Kp]); hold on
plot(massescon.Ms, [massescon.E; massescon.Kp],'--');
legend('ScaleE E', 'ScaleE Kp', 'ConsV E', 'ConsV Kp','Location','South');
xlabel('Total mass'); title('Energy and Kp');
subplot(223); plot(masses.Ms, [masses.totalwork' masses.comphases]); hold on;
plot(massescon.Ms, [massescon.totalwork' massescon.comphases], '--'); 
legend('total', 'CO', 'RB', 'PL', 'PO','Location','Best');
xlabel('Total mass'); ylabel('Work/step'); title('COM work per phase');
subplot(224); plot(masses.Ms, [masses.ds; masses.dsAbs], '-'); hold on;
plot(massescon.Ms, [massescon.ds; massescon.dsAbs], '--'); 
legend('ScaleE DS frac', 'ScaleE DS abs', 'ConsV DS frac', 'ConsV DS abs','Location','Best');
xlabel('Total mass'); ylabel('Double support'); title('Double support time');
% it seems that humans do swing the leg harder just as they increase their
% overall energy. 
%% ConsV Results: Min & max vertical forces, v^2 delta^2 predictor
% We find that vertical force excursion increases somewhat linearly
% with mass, and the M v^2 delta^2 predictor is pretty decent with mass.
clf;
subplot(221)
plot(masses.Ms, [masses.minmaxfy diff(masses.minmaxfy,1,2)])
hold on; plot(massescon.Ms, [massescon.minmaxfy diff(massescon.minmaxfy,1,2)], '--');
xlabel('Total mass'); ylabel('Force'); title('Min & max vert forces');
legend('min','max','diff','Location','East');
subplot(222); 
%plot(masses.Ms, masses.comydisp); hold on; % plot COM displacement
%plot(massescon.Ms, massescon.comydisp, '--'); set(gca,'ylim',[0 Inf]);
%xlabel('Total mass'); ylabel('Displacement'); title('Vert COM displ');
plot(0.5*masses.Ms .* masses.vminus.^2 .* masses.delta.^2, masses.totalwork); hold on
plot(0.5*massescon.Ms .* massescon.vminus.^2 .* massescon.delta.^2, massescon.totalwork, '--')
plot([0 0.2],[0 0.2],':'); axis equal; set(gca,'Xlim',[0 Inf],'YLim',[0 Inf]);
legend('ScaleE', 'ConsV','Location','Best');
xlabel('1/2 M v^2 delta^2'); ylabel('Total work');
subplot(223);
plot(masses.Ms, masses.delta*180/pi); set(gca,'ylim',[0 Inf]); hold on;
plot(massescon.Ms, massescon.delta*180/pi, '--');
xlabel('Total mass'); ylabel('Delta (deg)'); title('Delta vcom');
subplot(224);
plot(masses.Ms, masses.vminus); set(gca,'ylim',[0 Inf]); hold on;
plot(massescon.Ms, massescon.vminus, '--'); 
xlabel('Total mass'); ylabel('vminus'); title('vminus & predictor');
plot(masses.Ms, 0.5*masses.Ms .* masses.vminus.^2 .* masses.delta.^2,'r'); hold on
plot(masses.Ms, 0.5*massescon.Ms .* massescon.vminus.^2 .* massescon.delta.^2,'r--')
legend('ScaleE v-','ConsV v-', 'ScaleE Mv^2\delta^2','ConsV Mv^2\delta^2','Location','West')
end % if you want to see ConsV results

%% Examine walking at higher speeds and step lengths
speeds = 0.3:0.05:0.5;
steplens = 0.55:0.05:0.75;
w = wstar3;
clear ssgrid
i = 1; % fill out step lengths
for j = 1:length(steplens)
  [w,cnvrg] = findgaitspeedb(w, 'speed', speeds(i), 'steplength', steplens(j),'parmvary1','E','parmvary2','Kp');
  ssgrid.w(1,j) = w;
  ssgrid.cnvrg(1,j) = cnvrg;
end

j = 1; % fill out speeds
for i = 1:length(speeds)
  [w,cnvrg] = findgaitspeedb(w, 'speed', speeds(i), 'steplength', steplens(j),'parmvary1','E','parmvary2','Kp');
  ssgrid.w(i,j) = w;
  ssgrid.cnvrg(i,j) = cnvrg;
end

for i = 1:length(speeds)
  if ssgrid.cnvrg(i,1)
    w = ssgrid.w(i,1);
  else
    w = wstar3;
  end
  for j = 1:length(steplens)      
    [w,cnvrg] = findgaitspeedb(w, 'speed', speeds(i), 'steplength', steplens(j),'parmvary1','E','parmvary2','Kp');
    if ~cnvrg % didn't work, so let's go back to wstar3
      [w,cnvrg] = findgaitspeedb(wstar3, 'speed', speeds(i), 'steplength', steplens(j),'parmvary1','E','parmvary2','Kp');
      if ~cnvrg
        error('could not find gait');
      end
    end
    ssgrid.w(i,j) = w;
    ssgrid.cnvrg(i,j) = cnvrg;
    [grf,grt,cop,ssgrid.work = groundreactionforces( 
  end
end
ssgrid.speeds = speeds;
ssgrid.steplens = steplens;

clf; hold on;
for i = 1:length(speeds);
  for j = 1:length(steplens)
    if ssgrid.cnvrg(i,j)
      groundreactionforces(ssgrid.w(i,j));
    end
    pause
  end
end

clf; hold on;
for i = 1:length(speeds);
  for j = 1:length(steplens)
    if ssgrid.cnvrg(i,j)
      onestep(ssgrid.w(i,j));
    end
  end
end


%% Verify that bvp solver works
wb0 = bvpsearch(wstar3);
wb = findgaitb(wstar3, set(wstar3,'E',1.1,'Kp',12))
wb2 = findgaitb(wstar3, set(wstar3,'E',1.8,'Kp',12))

speedinfo = walkspeed(wstar3);
[w,cnvrg] = findgaitspeedb(wstar3, 'speed', speedinfo.speed, 'steplength', speedinfo.steplength,'parmvary1','E','parmvary2','Kp','info',2);

tic;[w,cnvrg,sol] = findgaitspeedb(wstar3, 'speed', 0.5, 'steplength', 0.7,'parmvary1','E','parmvary2','Kp','info',2);toc
tic;[w,cnvrg] = findgaitspeed(wstar3, 'speed', 0.5, 'steplength', 0.7,'parmvary1','E','parmvary2','Kp','info',2);toc







