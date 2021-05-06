
% Simulate cetain number of steps with certain Lgain conditions, 
% and then generate animation. 
%
% by Hansol Ryu, for Ryu&Kuo paper on central pattern generator. 
% tested on Matlab R2019a. 

addpath(fullfile(pwd,'..\..\src'))

load('result_Reference.mat');
% load conditions used for paper analysis 

maxTrial = 1;
noise_miss = [0 inf];

parms = get(w0,'parms');
g = parms.g; Mp = parms.Mp; M = parms.M;
L = parms.L; C = parms.C; R = parms.R; Ip = parms.Ip; Il = parms.Il;

Tst0 = get(w0, 'Tst');

seedInt = 59;
limit_max = 2000;
seed = randseed(seedInt, limit_max);
anim = 16;
est_time = 2; % time limit for one step
maxStep = 16;
plotRange = 1:maxStep;

ww = 0.005; % process noise deviation
vv = 0.10;   % sensory noise deviation
noiseVector = [MM_L\[ww;ww]; vv; vv];

G = eye(4);
Q_ = diag([0,0,noiseVector(1:2)']).^2;
R_ = diag(noiseVector(3:4)).^2;
N = zeros(4,2);

result_cell = cell(maxTrial, length(noise_miss)+1, 6);
% 1: num_steps, 2: step_len, 3: st_energy, 4:sw_energy, 5:st_time
% 6: isfall
result_num_fall = nan(maxTrial, length(noise_miss)+1,2); %1:num 2:time
result_num_steps = nan(maxTrial, length(noise_miss)+1);
result_dist = nan(maxTrial, length(noise_miss)+1);
result_time = nan(maxTrial, length(noise_miss)+1);
result_rmsErr = nan(maxTrial, length(noise_miss)+1);

result_stE = nan(maxTrial, length(noise_miss)+1);
result_swE = nan(maxTrial, length(noise_miss)+1);
result_stepVar = nan(maxTrial, length(noise_miss)+1);
result_motion = cell(maxTrial, length(noise_miss)+1);

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

% simulate for for maxStep-steps 
for trial_itr = 1:maxTrial
    % generate noise for a given trial % 1 8 9 12
    s = rng(seed(20), 'v5normal');
    noise_dt = est_time/32;
    noise_t = 0:noise_dt:est_time*100*3;
    
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
    
    %             noise_signal( round(1.5*est_time/noise_dt), 2) = 5;
    % %             for happy ending figure only!
    
    for noise_itr = 1:length(noise_miss)+1
        w = w0; % may not need...
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
        
        %                 L_norm(end) = L_norm(end-1)/5^3;
        
        w = set(w, 'Lgain', Lgain);
        
        t_global = 0;
        
        step_len = nan(maxStep, 1);
        est_err_squared = nan(maxStep, 3); %1:stance, 2:swing, 3:#of samples
        energy_stance = nan(maxStep, 2);
        energy_swing = nan(maxStep,2);
        stance_time = nan(maxStep,1);
        
        x0 = get(w, 'xstar');
        xe = [x0, ...
            x0+[noise_signal(end,3) noise_signal(end,4) 0 0]];%
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
                xi = xe; % reset, and add penalty
                w = set(w, 'internalState', 1);
                step_len(k) = target_sL;
            else
                step_len(k) = R*(xs(1,1)-xs(end,1))+...
                    (L-R)*(sin(xs(end,2))-sin(xs(end,1)));
            end
            t_global = t_global+te;
%             step_len(k) = R*(xs(1)+xi(1))+(L-R)*(sin(xe(1))+sin(xi(1)));
            
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
        result_motion{trial_itr, noise_itr} = motionData;
        
        result_num_steps(trial_itr, noise_itr) = k;
        result_dist(trial_itr, noise_itr) = sum(step_len(plotRange));
        result_time(trial_itr, noise_itr) = sum(stance_time(plotRange));
        result_stE(trial_itr, noise_itr) = sum(energy_stance(plotRange,1));
        result_swE(trial_itr, noise_itr) = sum(energy_swing(plotRange,1));
        result_stepVar(trial_itr, noise_itr) = std(step_len(plotRange));
        
        result_rmsErr(trial_itr, noise_itr) = ...
            ( sum(est_err_squared(:,1))/...
            sum(est_err_squared(:,3)) )^0.5;
        
        falls = result_cell{trial_itr, noise_itr,6};
        step_len(falls==1) = [];
        stance_time(falls==1) = [];
        energy_stance(falls==1,:) = [];
        energy_swing(falls==1,:) = [];
    end % end for:noise_itr
end

%% check walking performance for simulated steps 
L0 = 1;
L_norm(isinf(L_norm)) = max(L_norm(~isinf(L_norm)))*2;
L_norm(L_norm==0) = min(L_norm(L_norm>0))/2;

figure

subplot(4,1,1:2)
energy_mean = nanmean(result_stE,1) + ...
    coeff*nanmean(result_swE,1);
semilogx((L_norm(1:end-2)/L_norm(L0)), ...
    energy_mean(1:end-2), 'o-k', 'markerFaceColor', [0 0 1], 'markerEdge', 'none');
hold on
plot(0.6, energy_mean(end), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')
plot(2, energy_mean(end-1), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')

subplot(4,1,3)
stepVar_mean = nanmean(result_stepVar,1);
semilogx((L_norm(1:end-2)/L_norm(L0)), ...
    stepVar_mean(1:end-2), 'ok-', 'markerFaceColor', [1 0 0], 'markerEdge', 'none');
hold on
plot(0.6, stepVar_mean(end), 'ok', 'markerFaceColor', [1 0 0], 'markerEdge', 'none')
plot(2, stepVar_mean(end-1), 'ok', 'markerFaceColor', [1 0 0], 'markerEdge', 'none')
xlim([0.6 2])

ax=subplot(4,1,4);
numFalls = zeros(length(noise_miss)+1, maxTrial);
for noise_itr=1:length(noise_miss)+1
    for trial_itr=1:maxTrial
        numFalls(noise_itr, trial_itr) = ...
            sum(result_cell{trial_itr, noise_itr, 6}(plotRange));
    end
end
numFalls_mean = mean(numFalls, 2);
semilogx((L_norm(1:end-2)/L_norm(L0)), ...
    numFalls_mean(1:end-2), 'o-k', 'markerFaceColor', [0 0 1], 'markerEdge', 'none');
hold on
plot(0.6, numFalls_mean(end), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')
plot(2, numFalls_mean(end-1), 'ok', 'markerFaceColor', [0 0 1], 'markerEdge', 'none')

%% animate
walkedVel = result_dist./result_time;
Ltext = {'L^*_{lqe}', '\infty (FB)', '0 (FF)'};
for noise_itr = 1:3%1:5 %2%1:length(noise_miss)+1
    figure('units','normalized','outerposition',[0 0 1 1])
    
    motionData = result_motion{1, noise_itr}(plotRange,:);
    
    energy_stance = result_cell{1, noise_itr,3}(plotRange,:);
    energy_swing = result_cell{1, noise_itr,4}(plotRange,:);
    energy_total = energy_stance(:,1) + energy_swing(:,1)*coeff;
    
    falls = result_cell{1, noise_itr,6}(plotRange,:);
    
    animateOptimal(w0, 'drawFromData', motionData(1:end,1), ...
        'numsteps', plotRange(end)-plotRange(1)+1, 'falls', falls,'treadmill', false, ...
        'stepEnergy', energy_total, ...
        'modelLabel', ['L = ', Ltext{noise_itr}]);
end