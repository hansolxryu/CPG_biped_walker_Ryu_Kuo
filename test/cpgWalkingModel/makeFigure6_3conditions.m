
F1 = figure; hold on;
F2 = figure; hold on;
F3 = figure; hold on;
F4 = figure; hold on;
F5 = figure; hold on;
F6 = figure; hold on;

for cond = 1:3
    switch cond
        case 1
            load('result_Low.mat');
            col = [119 112 100]/255;
            mar = 'v';
        case 2
            load('result_Medium.mat');
            col = [206 0 54]/255;
            mar = 'o';
        case 3
            load('result_High.mat');
            col = [80 116 188]/255;
            mar = '^';
    end

COT_mean = nanmean(result_stE./result_dist) + ...
    coeff*nanmean(result_swE./result_dist);
COT_std = nanstd(result_stE./result_dist) + ...
    coeff*nanstd(result_swE./result_dist);

COT_nofall_mean = nanmean(nofall_result_stE./nofall_result_dist) + ...
    coeff*nanmean(nofall_result_swE./nofall_result_dist);
COT_nofall_std = nanstd(nofall_result_stE./nofall_result_dist) + ...
    coeff*nanstd(nofall_result_swE./nofall_result_dist);

% cm = colormap(parula(10));
% col = cm(cond*3,:);

figure(F1);
ax =gca; % subplot(12,1,2:6);
yyaxis left
hold on
set(ax, 'XScale', 'log')

plot((L_norm(1:end-2)/L_norm(L0)), COT_mean(1:end-2),...
    'marker',mar, 'markerFaceColor', col, 'markerEdge', 'none');
plot(2, COT_mean(end-1), ...
    'marker',mar, 'markerFaceColor', col, 'markerEdge', 'none')

% plot((L_norm(1:end-2)/L_norm(L0)), COT_nofall_mean(1:end-2),...
%     'marker',mar,  'markerFaceColor', col*0.8, 'markerEdge', 'none');
% plot([2 0.6], COT_nofall_mean(end-1:end), ...
%     'marker',mar,  'markerFaceColor', col*0.8, 'markerEdge', 'none')
ytickformat('%.2f')

ylabel('CoT')

semilogx([0.6;(L_norm(1:end-2)/L_norm(L0));2], ...
    ones(length(L_norm),1)*COT0, '-', 'color', col)
ylim([0.05 0.11])
yl=ylim;
yLength2to6 = ax.Position(4);
yyaxis right
semilogx((L_norm(1:5)/L_norm(L0)), ones(5,1)*100, '-','color', col)
ylim(yl./COT0*100)
yticks(100:20:200)
xlim([0.55 2.2])
xl = xlim;

figure(F2);
ax = gca;%subplot(12,1,1); % FF COT ends up being to high, on separate 
yLength1 = ax.Position(4);
yyaxis left
hold on
% plot(0.6, COT_mean(end), 'k')
semilogx(0.6, COT_mean(end), 'marker',mar, 'markerFaceColor', col, 'markerEdge', 'none')
text(0.7, COT_mean(end), ...
    [num2str(COT_mean(end),'%.2f'),'+-',num2str(COT_std(end),'%.2f'),...
    ',',num2str(COT_mean(end)*100/COT0,'%.1f'),'%']);
ylim([0.335, (yl(2)-yl(1))/yLength2to6*yLength1+0.335])
yl = ylim;
yticks(0.35:0.01:0.40)
yyaxis right
semilogx((L_norm(1:5)/L_norm(L0)), ones(5,1)*100*COT_mean(end)/COT0, '-k')
ylim(yl./COT0*100)
xlim(xl)
yticks([630, 640]);
set(gca, 'XScale', 'log')
yticks(600:20:800)

figure(F3);
% subplot(6,1,4)
stepVar_mean = nanmean(result_stepVar,1);
stepVar_std = nanstd(result_stepVar,1);

plot((L_norm(1:end-2)/L_norm(L0)), stepVar_mean(1:end-2), ...
    'marker',mar, 'markerFaceColor', col, 'markerEdge', 'none');
hold on
plot([2 0.6], stepVar_mean(end-1:end), ...
    'marker',mar, 'lineStyle','none', 'markerFaceColor', col, 'markerEdge', 'none')
ylim([0.04 0.08])
yticks(0.04:0.01:0.08)
ytickformat('%.2f')
xlim(xl)
ylabel('step length variation');
set(gca, 'XScale', 'log')

figure(F4);
% ax=subplot(6,1,5);
MTBF_mean = mean(MTBF_result, 2);
MTBF_std = std(MTBF_result');

plot((L_norm(1:end-2)/L_norm(L0)), MTBF_mean(1:end-2), ...
    'marker',mar, 'markerFaceColor', col, 'markerEdge', 'none');
hold on
plot([2 0.6], MTBF_mean(end-1:end), ...
    'marker',mar, 'lineStyle','none','markerFaceColor', col, 'markerEdge', 'none')

xlim(xl)
errorbar(1.6, 3, target_sL/target_v/2)
text(1.65,3, 'nominal step period')

yticks(0:2:12)
ylabel('mean time between falls');
set(gca, 'XScale', 'log')

figure(F5);
ax=gca;%subplot(12,1,11); % FF and FB has too high estimation error
RMS_mean = mean(result_rmsErr);
RMS_std = std(result_rmsErr);

hold on
plot([2 0.6], RMS_mean(end-1:end), ...    
    'marker',mar, 'lineStyle','none', 'markerFaceColor', col, 'markerEdge', 'none')
text(0.65, RMS_mean(end), [num2str(RMS_mean(end),'%.2f'),'+-',...
    num2str(RMS_std(end),'%.2f')]);
text(2.05, RMS_mean(end-1), num2str(RMS_mean(end-1),'%.2f'));
set(gca, 'XScale', 'log')

ylim([3 9])
xlim(xl)
ax.YAxis.MinorTick ='on';
yticks([2:2:8])
ax.YAxis.MinorTickValues = 2:9;
ytickformat('%.2f')

figure(F6);
ax=gca;%subplot(12,1,12);
hold on
plot((L_norm(1:end-2)/L_norm(L0)), RMS_mean(1:end-2), ...
    'marker', mar, 'markerFaceColor', col, 'markerEdge', 'none');xlim(xl)
set(gca, 'XScale', 'log')

xticks([0.6:0.2:1.6, 2])
ax.YAxis.MinorTick ='on';
yticks(0.08:0.02:0.15);
ytickformat('%.2f')
ax.YAxis.MinorTickValues = 0.07:0.01:0.15;
ylabel('estimation error');
xtickformat('%.1f')

end