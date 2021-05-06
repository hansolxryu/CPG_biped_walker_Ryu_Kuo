figure
load('result_Reference.mat')

COT_mean = nanmean(result_stE./result_dist) + ...
    coeff*nanmean(result_swE./result_dist);
COT_std = nanstd(result_stE./result_dist) + ...
    coeff*nanstd(result_swE./result_dist);

COT_nofall_mean = nanmean(nofall_result_stE./nofall_result_dist) + ...
    coeff*nanmean(nofall_result_swE./nofall_result_dist);
COT_nofall_std = nanstd(nofall_result_stE./nofall_result_dist) + ...
    coeff*nanstd(nofall_result_swE./nofall_result_dist);

% figure()
cm = colormap(parula(10));
col = cm(2,:);


ax = subplot(12,1,2:6);
yyaxis left
hold on
set(ax, 'XScale', 'log')

errorbar((L_norm(1:end-2)/L_norm(L0)), COT_mean(1:end-2),...
    COT_std(1:end-2), 'o-k', 'markerFaceColor', col, 'markerEdge', 'none');
errorbar(2, COT_mean(end-1), COT_std(end-1), ...
    'ok', 'markerFaceColor', col, 'markerEdge', 'none')

errorbar((L_norm(1:end-2)/L_norm(L0)), COT_nofall_mean(1:end-2),...
    COT_nofall_std(1:end-2), 'o-k', 'markerFaceColor', col*0.8, 'markerEdge', 'none');
errorbar([2 0.6], COT_nofall_mean(end-1:end), COT_nofall_std(end-1:end), ...
    'ok', 'markerFaceColor', col*0.8, 'markerEdge', 'none')
ytickformat('%.2f')

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
set(gca, 'color','none')


ax = subplot(12,1,1); % FF COT ends up being to high, on separate 
yLength1 = ax.Position(4);
yyaxis left
hold on
errorbar(0.6, COT_mean(end), COT_std(end), 'k')
semilogx(0.6, COT_mean(end), 'ok', 'markerFaceColor', col, 'markerEdge', 'none')
text(0.7, COT_mean(end), ...
    [num2str(COT_mean(end),'%.2f'),'+-',num2str(COT_std(end),'%.2f'),...
    ',',num2str(COT_mean(end)*100/COT0,'%.1f'),'%']);
ylim([0.335, (yl(2)-yl(1))/yLength2to6*yLength1+0.335])
yl = ylim;
yticks([0.33,0.34]);
yyaxis right
semilogx((L_norm(1:5)/L_norm(L0)), ones(5,1)*100*COT_mean(end)/COT0, '-k')
ylim(yl./COT0*100)
xlim(xl)
yticks([630, 640]);
ax = gca;
ax.YAxis(2).MinorTick ='on';
ax.YAxis(2).MinorTickValues = 625:5:635;
set(gca, 'XScale', 'log')
set(gca, 'color','none')

subplot(6,1,4)
stepVar_mean = nanmean(result_stepVar,1);
stepVar_std = nanstd(result_stepVar,1);

errorbar((L_norm(1:end-2)/L_norm(L0)), stepVar_mean(1:end-2), ...
    stepVar_std(1:end-2), 'ok-', 'markerFaceColor', col, 'markerEdge', 'none');
hold on
errorbar([2 0.6], stepVar_mean(end-1:end), ...
    stepVar_std(end-1:end),'ok', 'markerFaceColor', col, 'markerEdge', 'none')
ylim([0.04 0.068])
ytickformat('%.3f')
xlim(xl)
ylabel('step length variation');
set(gca, 'XScale', 'log')
set(gca, 'color','none')

ax=subplot(6,1,5);
MTBF_mean = mean(MTBF_result, 2);
MTBF_std = std(MTBF_result');

errorbar((L_norm(1:end-2)/L_norm(L0)), MTBF_mean(1:end-2), ...
    MTBF_std(1:end-2), 'o-k', 'markerFaceColor', col, 'markerEdge', 'none');
hold on
errorbar([2 0.6], MTBF_mean(end-1:end), ...
    MTBF_std(end-1:end),'ok', 'markerFaceColor', col, 'markerEdge', 'none')

xlim(xl)
errorbar(1.6, 3, target_sL/target_v/2)
text(1.65,3, 'nominal step period')

yticks(0:2:12)
ylabel('mean time between falls');
set(gca, 'XScale', 'log')
set(gca, 'color','none')

ax=subplot(12,1,11); % FF and FB has too high estimation error
RMS_mean = mean(result_rmsErr);
RMS_std = std(result_rmsErr);

hold on
errorbar([2 0.6], RMS_mean(end-1:end), RMS_std(end-1:end), ...    
    'ok', 'markerFaceColor', col, 'markerEdge', 'none')
text(0.65, RMS_mean(end), [num2str(RMS_mean(end),'%.2f'),'+-',...
    num2str(RMS_std(end),'%.2f')]);
text(2.05, RMS_mean(end-1), num2str(RMS_mean(end-1),'%.2f'));
set(gca, 'XScale', 'log')

ylim([2.5 4])
xlim(xl)
ax.YAxis.MinorTick ='on';
yticks([3 4])
ax.YAxis.MinorTickValues = 2.5:0.1:4;
ytickformat('%.2f')
set(gca, 'color','none')

ax=subplot(12,1,12);
hold on
errorbar((L_norm(1:end-2)/L_norm(L0)), RMS_mean(1:end-2), ...
    RMS_std(1:end-2), 'o-k', 'markerFaceColor', col, 'markerEdge', 'none');xlim(xl)
set(gca, 'XScale', 'log')

xticks([0.6:0.2:1.6, 2])
ax.YAxis.MinorTick ='on';
yticks(0.07:0.01:0.1);
ytickformat('%.2f')
ax.YAxis.MinorTickValues = 0.07:0.005:0.1;
ylabel('estimation error');
xtickformat('%.1f')
set(gca, 'color','none')