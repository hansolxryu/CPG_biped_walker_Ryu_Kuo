
% Get walking performances with various Lgains. 
% Simulate cetain number of steps with certain Lgain conditions, 
% and obtain walking performances in terms of energy, frequency of falls, 
% variability, estimation error. Repeat for several random noise profiles. 
% for paper, number of steps("maxStep") = 100; 
%            7 Lgain conditions(6 specified by "noise_miss", pure FF);
%            number of repeatition("maxTrial") = 20;
%
% by Hansol Ryu, for Ryu&Kuo paper on central pattern generator. 
% tested on Matlab R2019a. 

addpath(fullfile(pwd,'..\..\src'))

%% initialize %%
w = walksw2cpg;
Lgain = zeros(4,2);
w = set(w, 'Lgain', Lgain);

parms = get(w,'parms');
g = parms.g; Mp = parms.Mp; M = parms.M;
L = parms.L; C = parms.C; R = parms.R; Ip = parms.Ip; Il = parms.Il;

% initial guess %
w = set(w, 'Tst', -.03, 'Kp', .25, 'theta0', 0);
w = set(w, 'xstar', [0.33   -0.33   -.4  -.2]);

% get limit cycle condition %
[w, cnvrg] = gradsearch(w, [], 'info', 0, 'criterion', 1e-9);
if(~cnvrg), error('not a good initial w'); end

% further search for desired speed and step length %
w = findgaitspeed(w, ...
    'speed', 0.4, 'steplength', 0.55, ...
    'parmvary1', 'Tst', 'parmvary2', 'Kp');
w0 = w;
vec = stability(w); % eigen values should be < 1
if(max(abs(vec))>1), error('unstable gait'); end

% check limit cycle from findgaitspeed() %
x0 = get(w, 'xstar');
xe = [x0, x0];
[xe,te,xs,ts, ~] = onestep(w, xe(1:4), ...
    'x0_h', xe(5:8), 'anim', 50, 'ts_global', 0, 'internal', true);

target_sL = R*(xs(1,1)-xs(end,1))+...
                        (L-R)*(sin(xs(end,2))-sin(xs(end,1)));
target_v = target_sL/te;

w = set(w, 'theta0', xe(1));
w0 = w;

nominal_theta = [xs(:,1), xs(:,3)];
nominal_phi = [xs(:,2), xs(:,4)];

coeff = 1.00; %0.99;%.51;
%% check COT
% energy consumption = stance leg work + coeff * swing leg work
% check whether we get optimal CoT at a given step length

num_point = -11; % set negative value to skip this process
w_ = w0;

if(num_point>0)
    result_sw_E = nan(num_point, 1);
    result_st_E = nan(num_point, 1);
    result_dist = nan(num_point, 1);
    result_time = nan(num_point, 1);
    
    result_Tst =  nan(num_point, 1);
    result_Tsw =  nan(num_point, 1);
    
    portion = linspace(0.9,1.1,num_point);
    
    w_ = w0;
    for k=(num_point-1)/2:-1:1
        k
        w_ = findgaitspeed(w_, ...
            'speed', target_v, 'steplength', target_sL*portion(k), ...
            'parmvary1', 'Tst', 'parmvary2', 'Kp');
        
        [sL, Est, Esw, sT] = sim_n_steps(w_, 10);
        result_sw_E(k) = sum(Esw(:,1));
        result_st_E(k) = sum(Est(:,1));
        result_dist(k) = sum(sL);
        result_time(k) = sum(sT);
        result_Tst(k) =  get(w_, 'Tst');
        result_Tsw(k) =  get(w_, 'Kp');
    end
    w_ = w0;
    for k=(num_point+1)/2:1:num_point
        k
        w_ = findgaitspeed(w_, ...
            'speed', target_v, 'steplength', target_sL*portion(k), ...
            'parmvary1', 'Tst', 'parmvary2', 'Kp');
        
        [sL, Est, Esw, sT] = sim_n_steps(w_, 10);
        result_sw_E(k) = sum(Esw(:,1));
        result_st_E(k) = sum(Est(:,1));
        result_dist(k) = sum(sL);
        result_time(k) = sum(sT);
        result_Tst(k) =  get(w_, 'Tst');
        result_Tsw(k) =  get(w_, 'Kp');
    end
    
    draw_range = 1:num_point;
    
    figure
    subplot(2,1,1)
    plot(target_sL*portion(draw_range), result_sw_E(draw_range)./result_dist(draw_range), 'k:')
    hold on
    plot(target_sL*portion(draw_range), result_st_E(draw_range)./result_dist(draw_range), 'k-.')
    ylabel('CoT (dimensionless)')
    xlabel(['step length (walking spd = ', num2str(target_v), ')'])
    legend('swing Energy', 'stance Energy')
    
    yl = ylim();
    plot([target_sL, target_sL], [yl(1), yl(2)], 'k:')
    
    %     coeff = 0.12;
    subplot(2,1,2)
    
    hold off
    plot(target_sL*portion(draw_range), ...
        (result_sw_E(draw_range)*coeff + ...
        result_st_E(draw_range))./result_dist(draw_range), 'k.--');
    
    ylabel('CoT (dimensionless)')
    xlabel(['step length (walking spd = ', num2str(target_v), ')'])
end

%% setup lqe parameters, get optimal & syboptimal L gains %%

% linearization %
MM_L = zeros(2,2); rhs_L = zeros(2,2);
MM_L(1,1) = Il + M*(C*C) + M*(L*L) + Mp*(L*L);
MM_L(1,2) = C*L*M - M*(L*L);
MM_L(2,1) = MM_L(1,2);
MM_L(2,2) = Il - 2*C*L*M + M*(C*C) + M*(L*L);

rhs_L(1,:) = [-g*M*(-C + R)-g*M*(-L + R)-g*Mp*(-L + R), 0];
rhs_L(2,:) = [0, -g*(-C + L)*M];

A_0 = MM_L\rhs_L;

% numerically check linearization  %
ep = 0.000001;
parms2 = parms;
parms2.Kp = 0;
parms2.Tst = 0;
temp = fwalk(0, [ep 0 0 0 0 0 0 0]', w, parms2)/ep;
if(sum(temp(3:4) - A_0(:,1)>ep)), error('check A matrix'); end
temp = fwalk(0, [0 ep 0 0 0 0 0 0]', w, parms2)/ep;
if(sum(temp(3:4) - A_0(:,2)>ep)), error('check A matrix'); end

% lqe parameters %
A = zeros(4);
A(3:4, 1:2) = A_0;
A(1,3) = 1;
A(2,4) = 1;

ww = 0.005; % process noise deviation
vv = 0.10;   % sensory noise deviation

noiseVector = [MM_L\[ww;ww]; vv; vv];

G = eye(4);
Q_ = diag([0,0,noiseVector(1:2)']).^2;
% E{ww'} = Q, dx = Ax + Bu + Gw
R_ = diag(noiseVector(3:4)).^2;
% E{vv'} = R, y = Cx + Du + v
N = zeros(4,2);

C_matrix = [1 0 0 0; 0 1 0 0];
fprintf('number of unobservable modes is %d \n', ...
    length(A) - rank(obsv(A,C_matrix)));

% lqe parameters for suboptimal Lgain %
noise_miss = [-2,-0.5,0,0.25,0.4,inf];
L0 = 3; % index of optimum L, s.t. noise_miss=0
L_norm = zeros(length(noise_miss)+1, 1);

% get suboptimal Lgain %
for noise_itr = 1:length(noise_miss)+1
    if(noise_itr==length(noise_miss)+1||isnan(noise_miss(noise_itr)))
        Lgain = zeros(4,2); % FF
    elseif(isinf(noise_miss(noise_itr)))
        Lgain = inf;
    else
        % get gains from false noise information %
        [Lgain, P, E] = lqe(A, G, C_matrix, ...
            Q_*10^(noise_miss(noise_itr)), ...
            R_*10^(-noise_miss(noise_itr)), N);
    end
    L_norm(noise_itr) = norm(Lgain, 2);
end
L_norm'/L_norm(L0)
% -> large "noise_miss" = pretending to have large process noise -> more FB
%%

% save nominal parameters %
initial_walking = w0;
Tst0 = get(w0, 'Tst');
Kp0 = get(w0, 'Kp');

anim = 16;
est_time = 2; % time limit for one step
maxStep = 100;
maxTrial = 20;
seedInt = 59;
limit_max = 2000;
seed = randseed(seedInt, limit_max);

result_cell = cell(maxTrial, length(noise_miss)+1, 6);
% 1: num_steps, 2: step_len, 3: st_energy, 4:sw_energy, 5:st_time
% 6: isfall
result_num_fall = nan(maxTrial, length(noise_miss)+1,2); %1:num 2:time
result_num_steps = nan(maxTrial, length(noise_miss)+1);
result_dist = nan(maxTrial, length(noise_miss)+1);
result_time = nan(maxTrial, length(noise_miss)+1);

result_stE = nan(maxTrial, length(noise_miss)+1);
result_swE = nan(maxTrial, length(noise_miss)+1);
result_stepVar = nan(maxTrial, length(noise_miss)+1);

result_rmsErr = nan(maxTrial, length(noise_miss)+1);

nofall_result_dist = nan(maxTrial, length(noise_miss)+1);
nofall_result_time = nan(maxTrial, length(noise_miss)+1);
nofall_result_stE = nan(maxTrial, length(noise_miss)+1);
nofall_result_swE = nan(maxTrial, length(noise_miss)+1);

try    
    for trial_itr = 1:maxTrial
        % generate noise for a given trial %
        s = rng(seed(trial_itr), 'v5normal');
        noise_dt = est_time/32;
        noise_t = 0:noise_dt:est_time*maxStep*6;
        
        noise_ = normrnd(0,1,8*length(noise_t),1);
        noise_(abs(noise_)>3) = [];
        noise_signal= zeros(length(noise_t),4);
        
        for k=1:4 % cut out too bad noise
            noise_signal(:,k) = noise_(length(noise_t)*(k-1)+1:...
                length(noise_t)*k);
            noise_signal(:,k) = noise_signal(:,k)-...
                mean(noise_signal(:,k));
            noise_signal(:,k) = noise_signal(:,k)./...
                std(noise_signal(:,k))*noiseVector(k);
        end
        
        for noise_itr = 1:length(noise_miss)+1            
            w = initial_walking; % start with initial gait
            w = set(w, 'Tst', Tst0);
            w = set(w, 'Kp', Kp0);
            w = set(w, 'internalState', 1);
            
            curr_v = target_v;
            motionData = {}; % for animation
            
            [trial_itr, noise_itr]
            G = eye(4);
            Q_ = diag([0,0,noiseVector(1:2)']).^2;
            % E{ww'} = Q, dx = Ax + Bu + Gw
            R_ = diag(noiseVector(3:4)).^2;
            % E{vv'} = R, y = Cx + Du + v
            N = zeros(4,2); 
            % E{wv'} = N
            
            if(noise_itr==length(noise_miss)+1||isnan(noise_miss(noise_itr)))
                Lgain = zeros(4,2); % FF
            elseif(isinf(noise_miss(noise_itr)))
                Lgain = inf;
            else
                % get gains from false noise information %
                [Lgain, P, E] = lqe(A, G, C_matrix, ...
                    Q_*10^(noise_miss(noise_itr)), ...
                    R_*10^(-noise_miss(noise_itr)), N);
            end
            L_norm(noise_itr) = norm(Lgain, 2);
           
            
            w = set(w, 'Lgain', Lgain);
            
            t_global = 0;
            
            step_len = nan(maxStep, 1);
            est_err_squared = nan(maxStep, 3); %1:stance, 2:swing, 3:#of samples
            energy_stance = nan(maxStep, 2);
            energy_swing = nan(maxStep,2);
            stance_time = nan(maxStep,1);
            
            x0 = get(w, 'xstar');
            xe = [x0, ...
                x0+[noise_signal(end,3) noise_signal(end,4) 0 0]];
            % put some initial error to estimated angles
            
            xi = xe; % initial angle of each step. to get CoT
            %                 figure
            result_cell{trial_itr, noise_itr,6} = zeros(100,1); %fall?
            result_num_fall(trial_itr, noise_itr,1) = 0;
            
            for k=1:maxStep
                % send part of noise
                idx1 = max(1, ...
                    round(t_global/noise_dt-0.5*est_time/noise_dt));
                idx2 = min(length(noise_t), ...
                    round(t_global/noise_dt+1.5*est_time/noise_dt));
                
                noise_spline = spline(noise_t(idx1:idx2)', ...
                    noise_signal(idx1:idx2,:)');
                w = set(w, 'testParmMdl', ...
                    noise_spline);
                
                [xe,te,xs,ts, energies, match, int_match] = ...
                    onestep(w, xe(1:4), ...
                    'x0_h', xe(5:8), 'anim', anim, ...
                    'ts_global', t_global, 'internal', true);
                
                if(int_match(end)==-1)
                    w = set(w, 'internalState', -1);
                end
                
                motionData{k,1} = xs;
                motionData{k,2} = ts+t_global;
                
                if(isinf(Lgain))
                    noise_eval = ppval(noise_spline, t_global + ts);
                    noise_der=fnder(noise_spline,1);
                    noise_dot_eval=ppval(noise_der,t_global+ts);

                    xs(:,5:6) = xs(:,5:6) + noise_eval(3:4,:)';
                    xs(:,7:8) = xs(:,7:8) + noise_dot_eval(3:4,:)';
                end
                
                Torque_st_h = get(w, 'Tst')*ones(size(xs(:,7)));
                Torque_sw_h = -(get(w,'Kp')*(xs(:,6)));
                
                Torque_st1 = Torque_st_h;
                Torque_sw1 = Torque_sw_h;
                
                % torque might be applied to wrong leg
                Torque_st1(match==-1) = Torque_sw_h(match==-1);
                Torque_sw1(match==-1) = Torque_st_h(match==-1);
                
                motionData{k,3} = Torque_st_h;
                motionData{k,4} = Torque_sw_h;
                
                dW1 = (Torque_st1.*xs(:,3));
                dW2 = (Torque_sw1.*(xs(:,4)));
                
                Work_st_pos = cumtrapz(ts, dW1.*(dW1>0));
                Work_st_neg = cumtrapz(ts, dW1.*(dW1<0));
                
                Work_sw_pos = cumtrapz(ts, dW2.*(dW2>0));
                Work_sw_neg = cumtrapz(ts, dW2.*(dW2<0));
                
                if(isnan(xe))
                    error('something wrong, we should not reach here');
                end
                
                if(abs((pi/2-xe(1))*(xe(1)+pi/2))<1e-3 || ...
                        abs((pi/2-xe(5))*(xe(5)+pi/2))<1e-3 )
                    % fall! %
                    result_num_fall(trial_itr, noise_itr,1) = ...
                        result_num_fall(trial_itr, noise_itr,1)+1;
                    result_num_fall(trial_itr, noise_itr,2) = ...
                        t_global+te;
                    
                    result_cell{trial_itr, noise_itr,6}(k) = 1;
                    
                    % restart
                    xe = [x0, ...
                        x0+[noise_signal(end-k,3) noise_signal(end-k,4) 0 0]];
                    % give some initial error to estimated angles
                    xi = xe; % re-set, may add penalty here if desired
                    w = set(w, 'internalState', 1);    
                    
                    step_len(k) = target_sL;
                else
                    step_len(k) = R*(xs(1,1)-xs(end,1))+...
                        (L-R)*(sin(xs(end,2))-sin(xs(end,1)));
                end
                
                t_global = t_global+te;
                
                est_err_squared(k,1) = sum((xs(:,5)-xs(:,1)).^2)+...
                    sum((xs(:,6)-xs(:,2)).^2)+...
                    sum((xs(:,7)-xs(:,3)).^2)+...
                    sum((xs(:,8)-xs(:,4)).^2); 
                est_err_squared(k,3) = size(xs,1);
                
                energy_stance(k,:) = [Work_st_pos(end), Work_st_neg(end)];
                energy_swing(k,:) = [Work_sw_pos(end), Work_sw_neg(end)];
                stance_time(k) = te;
                
                curr_v = step_len(k)/te;
                w = set(w, 'Tst', ...
                    Tst0*(1-(sum(step_len(1:k))-target_v*t_global)*0.1));
                xi = xe;
            end % end for: step_itr
            
            result_cell{trial_itr, noise_itr,1} = k;
            result_cell{trial_itr, noise_itr,2} = step_len;
            result_cell{trial_itr, noise_itr,3} = energy_stance;
            result_cell{trial_itr, noise_itr,4} = energy_swing;
            result_cell{trial_itr, noise_itr,5} = stance_time;
            
            result_num_steps(trial_itr, noise_itr) = k;
            result_dist(trial_itr, noise_itr) = sum(step_len);
            result_time(trial_itr, noise_itr) = sum(stance_time);
            result_stE(trial_itr, noise_itr) = sum(energy_stance(:,1));
            result_swE(trial_itr, noise_itr) = sum(energy_swing(:,1));
            result_stepVar(trial_itr, noise_itr) = std(step_len);
            
            result_rmsErr(trial_itr, noise_itr) = ...
                ( sum(est_err_squared(:,1))/...
                sum(est_err_squared(:,3)) )^0.5;
            
            falls = result_cell{trial_itr, noise_itr,6};
            step_len(falls==1) = [];
            stance_time(falls==1) = [];
            energy_stance(falls==1,:) = [];
            energy_swing(falls==1,:) = [];
            
            nofall_result_dist(trial_itr, noise_itr) = sum(step_len);
            nofall_result_time(trial_itr, noise_itr) = sum(stance_time);
            nofall_result_stE(trial_itr, noise_itr) = sum(energy_stance(:,1));
            nofall_result_swE(trial_itr, noise_itr) = sum(energy_swing(:,1));
        end % end for:noise_itr
    end
catch ME
    display(ME);
    fprintf('%.3f %.3f\n', ww, vv);
end


%% calculate COT steady gait under no noise, at the outcome walking speeds

% to make sure optimal COT is not assosiated with slower walking speed
% In fact, optimal gait yeilds fastest walking speed.
mean_walk_spd = mean(result_dist./result_time);

% calculate nominal COT at each outcome walking speed.
% some FF conditions may not converge because it's walking at too low speed
stCOT_0noise = nan(maxTrial, length(noise_miss)+1);
swCOT_0noise = nan(maxTrial, length(noise_miss)+1);
calc_noFall_COT = false; % this takes time, turn off if not needed
if(calc_noFall_COT) 
    for trial_itr = 1:maxTrial
        for noise_itr = 1:length(noise_miss)+1
            spd_ = result_dist(trial_itr, noise_itr)./...
                result_time(trial_itr, noise_itr);
            if(isnan(spd_)), break; end;
            sl_ = result_dist(trial_itr, noise_itr)/...
                result_num_steps(trial_itr, noise_itr);
            
            w_ = findgaitspeed(w0, ...
                'speed', spd_, 'steplength', sl_, ...
                'parmvary1', 'Tst', 'parmvary2', 'Kp');
            %                 w_ = w0;n
            
            [sL, Est, Esw, sT] = sim_n_steps(w_, 10);
            stCOT_0noise(trial_itr, noise_itr) = sum(Esw(:,1))/sum(sL);
            swCOT_0noise(trial_itr, noise_itr) = sum(Est(:,1))/sum(sL);
        end
    end
    figure
    subplot(2,1,1)
    semilogx(L_norm/L_norm(L0), mean_walk_spd, '.', 'markerSize', 10);
    ylabel('result walking speed')
    subplot(2,1,2)    
    semilogx(L_norm/L_norm(L0), ...
        mean(stCOT_0noise)+mean(swCOT_0noise)*coeff, '.', 'markerSize', 10)
    ylabel('associated COT at the outcome speed')
    xlabel('L gain')    
    % better estimation -> higher walking speed -> higher nominal COT.
    % but COT was lowest at optimal Lgain.
    % -> nominal COT at outcome walkings speed doesn't explain COT trend.
    % -> supporting "better estimation yeilds better efficiency". 
end

%% plot result %%

[sL, Est, Esw, sT] = sim_n_steps(w0, 10);
COT0 = (sum(Est(:,1))+sum(Esw(:,1))*coeff)/sum(sL);

L_norm(isinf(L_norm)) = max(L_norm(~isinf(L_norm)))*2;
L_norm(L_norm==0) = min(L_norm(L_norm>0))/2;

COT_mean = nanmean(result_stE./result_dist) + ...
    coeff*nanmean(result_swE./result_dist);
COT_nofall_mean = nanmean(nofall_result_stE./nofall_result_dist) + ...
    coeff*nanmean(nofall_result_swE./nofall_result_dist);

figure()
ax = subplot(12,1,2:6);
yyaxis left
semilogx((L_norm(1:end-2)/L_norm(L0)), ...
    COT_mean(1:end-2), 'o-k', 'markerFaceColor', [0 0 1], 'markerEdge', 'none');
hold on
plot(2, COT_mean(end-1), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')

semilogx((L_norm(1:end-2)/L_norm(L0)), ...
    COT_nofall_mean(1:end-2), 'o-k', 'markerFaceColor', [1 0 0], 'markerEdge', 'none');
plot(2, COT_nofall_mean(end-1), 'ok', 'markerFaceColor', [1 0 0], 'markerEdge', 'none')
plot(0.6, COT_nofall_mean(end), 'ok', 'markerFaceColor', [1 0 0], 'markerEdge', 'none')
ylabel('CoT')

semilogx([0.6;(L_norm(1:end-2)/L_norm(L0));2], ...
    ones(length(L_norm),1)*COT0, '-k')
ylim([0.05 0.11])
yl=ylim;
yLength2to6 = ax.Position(4);
yyaxis right
semilogx((L_norm(1:5)/L_norm(L0)), ones(5,1)*100, '-k')
ylim(yl./COT0*100)

xlim([0.55 2.2])
xl = xlim;

ax = subplot(12,1,1); % FF COT ends up being to high, on separate 
yLength1 = ax.Position(4);
yyaxis left
semilogx(0.6, COT_mean(end), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')
hold on
text(0.7, COT_mean(end), ...
    [num2str(COT_mean(end)),',',num2str(COT_mean(end)*100/COT0),'%']);
ylim([0.335, (yl(2)-yl(1))/yLength2to6*yLength1+0.335])
yl = ylim;
yticks([0.33,0.34]);
yyaxis right
semilogx((L_norm(1:5)/L_norm(L0)), ones(5,1)*100*COT_mean(end)/COT0, '-k')
ylim(yl./COT0*100)
xlim(xl)
yticks(630);
ax = gca;
ax.YAxis(2).MinorTick ='on';
ax.YAxis(2).MinorTickValues = 625:5:635;

subplot(6,1,4)
stepVar_mean = nanmean(result_stepVar,1);
semilogx((L_norm(1:end-2)/L_norm(L0)), ...
    stepVar_mean(1:end-2), 'ok-', 'markerFaceColor', [0 1 1], 'markerEdge', 'none');
hold on
plot(0.6, stepVar_mean(end), 'ok', 'markerFaceColor', [0 1 1], 'markerEdge', 'none')
plot(2, stepVar_mean(end-1), 'ok', 'markerFaceColor', [0 1 1], 'markerEdge', 'none')
ylim([0.044 0.068])
ytickformat('%.3f')
xlim(xl)
ylabel('step length variation');

ax=subplot(6,1,5);
MTBF_result = zeros(length(noise_miss)+1, maxTrial);
for noise_itr=1:length(noise_miss)+1%5
    for trial_itr=1:maxTrial
        temp = result_cell{trial_itr, noise_itr, 5};
        tempFall = result_cell{trial_itr, noise_itr, 6};
        temp(tempFall==1)=nan;
        
        temp2 = zeros(sum(tempFall),1);
        curr = 1;
        
        for step_itr=1:maxStep
            if(curr > sum(tempFall)),
                break;
            elseif tempFall(step_itr)==1 % reached fall
                curr = curr+1;
            else
                temp2(curr) = temp2(curr)+temp(step_itr);
            end
        end
        MTBF_result(noise_itr, trial_itr) = mean(temp2);
    end
end
MTBF_mean = mean(MTBF_result, 2);
% mean(result_num_fall(:,:,1))
semilogx((L_norm(1:end-2)/L_norm(L0)), ...
    MTBF_mean(1:end-2), 'o-k', 'markerFaceColor', [0 0 1], 'markerEdge', 'none');
hold on
plot(0.6, MTBF_mean(end), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')
plot(2, MTBF_mean(end-1), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')
xlim(xl)
errorbar(1.6, 3, target_sL/target_v/2)
yticks(0:2:12)
ylabel('number of falls');

ax=subplot(12,1,11); % FF and FB has too high estimation error
RMS_mean = mean(result_rmsErr);
hold off
semilogx(0.6, RMS_mean(end), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')
text(0.6, RMS_mean(end)-0.5, num2str(RMS_mean(end)));
hold on
semilogx(2, RMS_mean(end-1), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')
text(2, RMS_mean(end-1)-0.5, num2str(RMS_mean(end-1)));
ylim([3 4])
xlim(xl)
ax.YAxis.MinorTick ='on';
ylim([2.5 4])
yticks(3:4)
ax.YAxis.MinorTickValues = 2.5:0.1:4;

ax=subplot(12,1,12);
hold off
semilogx((L_norm(1:end-2)/L_norm(L0)), ...
    RMS_mean(1:end-2), 'o-k', 'markerFaceColor', [0 0 1], 'markerEdge', 'none');xlim(xl)
xticks([0.6:0.2:1.6, 2])

ax.YAxis.MinorTick ='on';
ax.YAxis.MinorTickValues = 0.07:0.005:0.1;
ylabel('estimation error');


MSBF_mean = zeros(length(noise_miss)+1, maxTrial);
% mean step numbers btw falls
for noise_itr=1:length(noise_miss)+1%5
    for trial_itr=1:maxTrial
    MSBF_mean(noise_itr, trial_itr) = ...
        find(result_cell{trial_itr,noise_itr,6}==1, 1, 'last')/...
        sum(result_cell{trial_itr,noise_itr,6}==1) - ...
        1;
    end
end
MSBF_mean = mean(MSBF_mean,2)';


numFall_mean = zeros(length(noise_miss)+1, maxTrial);
% mean step numbers btw falls
for noise_itr=1:length(noise_miss)+1%5
    for trial_itr=1:maxTrial
    numFall_mean(noise_itr, trial_itr) = ...
        sum(result_cell{trial_itr,noise_itr,6}==1);
    end
end
numFall_mean = mean(numFall_mean,2)';


walk_spd_mean = nanmean(result_dist./result_time,1);
walk_sl_mean = nanmean(result_dist./result_num_steps,1);
zeroNoiseCOT_mean = mean(swCOT_0noise*coeff+stCOT_0noise);

% summary measures
bigResTable = [COT_mean; COT_nofall_mean; stepVar_mean; ...
    MTBF_mean'; MSBF_mean;
    walk_spd_mean; walk_sl_mean; zeroNoiseCOT_mean];
bigResTable(:,[end 1:5 end-1])


%%
function [step_len, energy_stance, energy_swing, ...
    stance_time] = sim_n_steps(w,n)
t_global = 0;

step_len = nan(n, 1);
energy_stance = nan(n, 2);
energy_swing = nan(n,2);
stance_time = nan(n,1);

L = get(w, 'L');
R = get(w, 'R');
x0 = get(w, 'xstar');
anim = 50;

xe = [x0, x0];%
xi = xe;

left = [];
right = [];
for k=1:n
    
    [xe,te,xs,ts, energies, match] = onestep(w, xe(1:4), ...
        'x0_h', xe(1:4), 'anim', anim, 'ts_global', t_global, 'internal', false);
    
    xs = [xs, xs];
    Torque_st_h = get(w, 'Tst')*ones(size(xs(:,7)));
    Torque_sw_h = -(get(w,'Kp')*(xs(:,6)));
    
    Torque_st1 = Torque_st_h;
    Torque_sw1 = Torque_sw_h;
    
    % torque might be applied to wrong leg (in FF case), so correct that 
    Torque_st1(match==-1) = Torque_sw_h(match==-1);
    Torque_sw1(match==-1) = Torque_st_h(match==-1);
    
    dW1 = (Torque_st1.*xs(:,3));
    dW2 = (Torque_sw1.*(xs(:,4)));
    
    Work_st_pos = cumtrapz(ts, dW1.*(dW1>0));
    Work_st_neg = cumtrapz(ts, dW1.*(dW1<0));
    
    Work_sw_pos = cumtrapz(ts, dW2.*(dW2>0));
    Work_sw_neg = cumtrapz(ts, dW2.*(dW2<0));
    
    if(isnan(xe)) % fall over
        break;
    end
    t_global = t_global+te;
    
    step_len(k) = R*(xs(1,1)-xs(end,1))+...
        (L-R)*(sin(xs(end,2))-sin(xs(end,1)));

    energy_stance(k,:) = [Work_st_pos(end), Work_st_neg(end)];
    energy_swing(k,:) = [Work_sw_pos(end), Work_sw_neg(end)];
    stance_time(k) = te;
    
    xi = xe;
end

end