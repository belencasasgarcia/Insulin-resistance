function [] = plotExperimentalData_lowHC(EXPDATA)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

color =  [0, 0.4470, 0.7410];         % Blue
marker='o';
markersize=6;
fontsize=14;

figure()

title ('Experimental data')

subplot(3,2,1)
errorbar(EXPDATA.time{7},EXPDATA.mean{7},EXPDATA.SD{7},'Color', color,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Glucose concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);
box off

title('Liver + islets, 11 mM','FontSize',16)

subplot(3,2,2)
errorbar(EXPDATA.time{8},EXPDATA.mean{8},EXPDATA.SD{8},'Color', color,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Insulin concentration (\muU/mL)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);
box off

title('Liver + islets, 11 mM','FontSize',16)

% subplot(3,2,5)
% errorbar(EXPDATA.time{5},EXPDATA.mean{5},EXPDATA.SD{5},'Color', color,...
%     'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)
% 
% xlabel ('Time (h)','FontSize',fontsize)
% ylabel ('Glucose concentration (mM)','FontSize',fontsize)
% set(gca,'TickDir','out','FontSize',fontsize);
% box off
% 
% title('Liver + islets, 5.5 mM','FontSize',16)
% 
% 
% subplot(3,2,6)
% errorbar(EXPDATA.time{6},EXPDATA.mean{6},EXPDATA.SD{6},'Color', color,...
%     'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)
% 
% xlabel ('Time (h)','FontSize',fontsize)
% ylabel ('Insulin concentration (\muU/mL)','FontSize',fontsize)
% set(gca,'TickDir','out','FontSize',fontsize);
% box off
% 
% title('Liver + islets, 5.5 mM','FontSize',16)

end

