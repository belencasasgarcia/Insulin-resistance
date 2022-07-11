function [] = plotExperimentalData(EXPDATA)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

color =  [0, 0.4470, 0.7410];         % Blue
color11mM =  [0, 0.4470, 0.7410];         % Blue
color5p5mM = [0.6350, 0.0780, 0.1840];    % Red

marker='o';
markersize=6;
fontsize=14;

figure()

title ('Experimental data')

subplot(3,2,1)
errorbar(EXPDATA.time{1},EXPDATA.mean{1},EXPDATA.SD{1},'Color', color,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Glucose concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);
box off

title('Liver + islets, 11 mM','FontSize',16)

subplot(3,2,2)
errorbar(EXPDATA.time{2},EXPDATA.mean{2},EXPDATA.SD{2},'Color', color,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Insulin concentration (\muU/mL)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);
box off

title('Liver + islets, 11 mM','FontSize',16)

subplot(3,2,3)
errorbar(EXPDATA.time{3},EXPDATA.mean{3},EXPDATA.SD{3},'Color', color,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Glucose concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);
box off

title('Only liver','FontSize',16)

subplot(3,2,4)
errorbar(EXPDATA.time{4},EXPDATA.mean{4},EXPDATA.SD{4},'Color', color,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Insulin concentration (\muU/mL)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);
box off

title('Only liver','FontSize',16)

subplot(3,2,5)
errorbar(EXPDATA.time{5},EXPDATA.mean{5},EXPDATA.SD{5},'Color', color,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Glucose concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);
box off

title('Liver + islets, 5.5 mM','FontSize',16)


subplot(3,2,6)
errorbar(EXPDATA.time{6},EXPDATA.mean{6},EXPDATA.SD{6},'Color', color,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Insulin concentration (\muU/mL)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);
box off

title('Liver + islets, 5.5 mM','FontSize',16)

%% Validation data

% Glucose data 

figure()

title ('Validation data, Glucose')

subplot(2,2,1)
errorbar(EXPDATA.time{7},EXPDATA.mean{7},EXPDATA.SD{7},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Glucose concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);

xlim([288 360])

box off

title('11 mM, dose A','FontSize',16)

subplot(2,2,2)
errorbar(EXPDATA.time{8},EXPDATA.mean{8},EXPDATA.SD{8},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Glucose concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);

xlim([288 360])

box off

title('5.5 mM, dose A','FontSize',16)

subplot(2,2,3)
errorbar(EXPDATA.time{9},EXPDATA.mean{9},EXPDATA.SD{9},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Glucose concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);

xlim([288 360])

box off

title('11 mM, dose B','FontSize',16)

subplot(2,2,4)
errorbar(EXPDATA.time{10},EXPDATA.mean{10},EXPDATA.SD{10},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Glucose concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);

xlim([288 360])

box off

title('5.5 mM, dose B','FontSize',16)

%% Insulin data

figure()

title ('Validation data, Insulin')

subplot(2,2,1)
errorbar(EXPDATA.time{11},EXPDATA.mean{11},EXPDATA.SD{11},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Insulin concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);

xlim([288 360])

box off

title('11 mM, dose A','FontSize',16)

subplot(2,2,2)
errorbar(EXPDATA.time{12},EXPDATA.mean{12},EXPDATA.SD{12},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Insulin concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);

xlim([288 360])

box off

title('5.5 mM, dose A','FontSize',16)

subplot(2,2,3)
errorbar(EXPDATA.time{13},EXPDATA.mean{13},EXPDATA.SD{13},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Insulin concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);

xlim([288 360])

box off

title('11 mM, dose B','FontSize',16)

subplot(2,2,4)
errorbar(EXPDATA.time{14},EXPDATA.mean{14},EXPDATA.SD{14},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',markersize)

xlabel ('Time (h)','FontSize',fontsize)
ylabel ('Insulin concentration (mM)','FontSize',fontsize)
set(gca,'TickDir','out','FontSize',fontsize);

xlim([288 360])

box off

title('5.5 mM, dose B','FontSize',16)

end

