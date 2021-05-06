%% Load carriage study
% Using walksw2dhlc, the double-humped simplest walker
% with load carriage.
% Walksw2dhlc contains a parameter M for body mass, and computes all forces
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

wstar = walksw2dhlc('normal'); % this is derived from walksw2 and walks the same

% demonstrate finding a gait with a slightly larger total mass
w = wstar; Estar = get(w, 'E');
%w2 = findgait(w, set(w,'M', 1.01)) 

% Here's an alternate gait with slightly higher energy
% The general trends turn out the same as the nominal gait above
%w=set(walksw2dhlc, 'xstar', [ -0.378013 0.9258 -0.245515 0.466164 -0.112335 0 ], ...
%  'Kp', 7.417, 'Kleg', 30, 'steplength', 0.620297, 'E', 1.07, 'M', 1);
%wstar = gradsearch(w); Estar = get(wstar, 'E');

walksw2dhlcDirectory = fileparts(which('walksw2dhlc'));
%% Parameter study increasing mass, scaling E with mass
% This turns out to cause step length to increase slightly, but
% overall speed stays about the same because step frequency goes down a
% bit. This is not too different from what happens in experiment.
Ms = 1:0.05:1.4; w = wstar;
if exist([walksw2dhlcDirectory '/masses.mat']) % pre-saved results
  load([walksw2dhlcDirectory '/masses.mat'])
else % do the study
clear masses
for i=1:length(Ms)
  [w,cnvrg] = findgait(w, set(w, 'M', Ms(i), 'E', Estar*Ms(i)));
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
  toTime = find(ies == 2, 1);     % toe-off time
  masses.ds(i) = ts(toTime) / tc; % double support as fraction of stride
  masses.dsAbs(i) = ts(toTime);   % double support as an absolute
  [grf,grt,cop,comworkrate,finfo] = groundreactionforces(ts, xs, masses.w(i));
  masses.comphases(i,:) = finfo.comphases; % CO, RB, PL, PO in array
  masses.totalwork(i) = finfo.totalwork;   % total negative work
  masses.dswork(i) = finfo.dswork;         % double support work 
  masses.minmaxfy(i,:) = finfo.minmax';    % includes min and max vert forces
  masses.comydisp(i) = max(xs(:,2))-min(xs(:,2)); % vert COM displacement
  masses.delta(i) = finfo.delta;
  masses.vminus(i) = finfo.vminus;
%plot(ts, comworkrate);
end
plot(masses.Ms, masses.comphases) % show that work per phase increases with mass
xlabel('Total mass'); ylabel('COM work'); legend('CO','RB','PL','PO');
% save masses masses
end % do the study
%% Parameter study keeping speed and step length unchanged
% Now do a study varying mass while keeping step length and frequency the
% same. This is found by varying total energy E, and hip spring Kp to
% maintain the nominal gait. This requires a slightly increasing Kp 
if exist([walksw2dhlcDirectory '/massescon.mat']) % pre-saved results
  load([walksw2dhlcDirectory '/massescon.mat']);
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
  toTime = find(ies == 2, 1);        % toe-off time
  massescon.ds(i) = ts(toTime) / tc; % double support as fraction of stride
  massescon.dsAbs(i) = ts(toTime);   % double support as an absolute
  [grf,grt,cop,comworkrate,finfo] = groundreactionforces(ts, xs, masses.w(i));
  massescon.comphases(i,:) = finfo.comphases; % CO, RB, PL, PO in array
  massescon.totalwork(i) = finfo.totalwork;   % total negative work
  massescon.dswork(i) = finfo.dswork;         % double support work 
  massescon.minmaxfy(i,:) = finfo.minmax';    % includes min and max vert forces
  massescon.comydisp(i) = max(xs(:,2))-min(xs(:,2)); % vert COM displacement
  massescon.delta(i) = finfo.delta;
  massescon.vminus(i) = finfo.vminus;
%plot(ts, comworkrate);
end
%plot(massescon.Ms, massescon.comphases) % show work increasing per phase
%plot(massescon.Ms, [massescon.ds; massescon.dsAbs]) % double support increases with mass
% save massescon massescon
end % do the study
%% Results: Speed & step length, energy, Kp, COM work per phase,
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

%% Results: COM motion, ground reaction forces, COM work rate, and
% minimum and maximum vertical forces.
% Increasing mass causes more COM motion, more excursion of vertical
% forces, more peaks in COM work rate, and a taller hodograph
% without much change in horizontal COM velocity.
clf; 
for i = 1:length(Ms)
  subplot(221); hold on; xlabel('COM x');  ylabel('COM y'); title('COM motion');
  [xc,tc,xs,ts] = onestep(masses.w(i)); plot(xs(1:end-1,1), xs(1:end-1,2));
  subplot(222); hold on; xlabel('time'); ylabel('Forces'); title('Forces');
  [grf,grt,cop,comworkrate,finfo] = groundreactionforces(ts, xs, masses.w(i));
  plot(ts, grf);   
  subplot(223); hold on; xlabel('Total mass'); ylabel('COM work rate'); title('COM work rate');
  plot(ts, comworkrate); 
  subplot(224); hold on; plot(xs(:,4), xs(:,5)); xlabel('COM x vel'); ylabel('COm y vel'); title('COM Hodograph');
  axis equal
end

%% Results: Min & max vertical forces, v^2 delta^2 predictor
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

