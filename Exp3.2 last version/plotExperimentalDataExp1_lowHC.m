% Define colors

global EXPDATA

colorhyper_highHC = [204, 0, 51]/255;                % Hyperglycemia, 50 µM cortisone
colorhyper_lowHC = [0, 64, 128]/255;                 % Hyperglycemia, 10 nM cortisone
colornormo_lowHC = [149, 214, 246]/255;              % Normoglycemia, 10 nM cortisone
colornormo_highHC = [235, 174, 9]/255;	             % Normoglycemia, 50 µM cortisone
marker='o';    

% Plot glucose data
figure()

% Glucose, GTT d1-3

subplot(2,2,1)

% Hyperglycemia, high HC
errorbar(EXPDATA.time{1}(1:4,1), EXPDATA.mean{1}(1:4,1), EXPDATA.SD{1}(1:4,1),'Color', colorhyper_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_highHC, 'MarkerSize',8)
hold on

%Low HC

% Hyperglycemia, low HC
errorbar(EXPDATA.time{7}(1:4,1), EXPDATA.mean{7}(1:4,1), EXPDATA.SD{7}(1:4,1),'Color', colorhyper_lowHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_lowHC, 'MarkerSize',8)
hold on

% Normoglycemia, low HC, with or /wo GTT
errorbar(EXPDATA.time{9}(1:4,1), EXPDATA.mean{9}(1:4,1), EXPDATA.SD{9}(1:4,1),'Color', colornormo_lowHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_lowHC, 'MarkerSize',8)
hold on

xticks([0 8 24 48 288 296 312 336 360])
ylim([0 14])

xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

subplot(2,2,2)

% Hyperglycemia, high HC
errorbar(EXPDATA.time{1}(1:end,1), EXPDATA.mean{1}(1:end,1), EXPDATA.SD{1}(1:end,1),'Color', colorhyper_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_highHC, 'MarkerSize',8)
hold on

% % Normoglycemia, high HC, with GTT
% errorbar(EXPDATA.time{3}(5:end,1), EXPDATA.mean{3}(5:end,1), EXPDATA.SD{3}(5:end,1),'Color', colornormo_highHC,...
%     'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_highHC, 'MarkerSize',8)
% hold on

% Normoglycemia, high GC, with GTT

% Normoglycemia, high HC with GTT
errorbar(EXPDATA.time{5}(:,1), EXPDATA.mean{5}(:,1), EXPDATA.SD{5}(:,1),'Color', colornormo_highHC,...
     'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_highHC, 'MarkerSize',8)
hold on


%Low HC

% Hyperglycemia, low HC
errorbar(EXPDATA.time{7}(1:end,1), EXPDATA.mean{7}(1:end,1), EXPDATA.SD{7}(1:end,1),'Color', colorhyper_lowHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_lowHC, 'MarkerSize',8)
hold on

%  Normoglycemia, low HC, with GTT
errorbar(EXPDATA.time{9}(1:end,1), EXPDATA.mean{9}(1:end,1), EXPDATA.SD{9}(1:end,1),'Color', colornormo_lowHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_lowHC, 'MarkerSize',8)
hold on

% % Normoglycemia, low HC, w/o GTT
% errorbar(EXPDATA.time{11}(5:end,1), EXPDATA.mean{11}(5:end,1), EXPDATA.SD{11}(5:end,1),'Color', colornormo_lowHC,...
%     'LineStyle', '--', 'Marker',marker,'MarkerFaceColor',colornormo_lowHC, 'MarkerSize',8)
% hold on
% ylim([0 14])

xticks([0 8 24 48 288 296 312 336 360])

xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

legend('Hyperglycemia, High cortisone', 'Normoglycemia, High cortisone, with GTT',...
    'Hyperglycemia, Low cortisone',...
    'Normoglycemia, Low cortisone, with GTT')

    
%% Plot insulin data
% Hyperglycemia, high HC

%Insulin, GTT d1-3

subplot(2,2,3)
errorbar(EXPDATA.time{2}(1:4,1), EXPDATA.mean{2}(1:4,1), EXPDATA.SD{2}(1:4,1),'Color', colorhyper_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_highHC, 'MarkerSize',8)
hold on


%Low HC

% Hyperglycemia, low HC
errorbar(EXPDATA.time{8}(1:4,1), EXPDATA.mean{8}(1:4,1), EXPDATA.SD{8}(1:4,1),'Color', colorhyper_lowHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_lowHC, 'MarkerSize',8)
hold on

% Normoglycemia, low HC, with or /wo GTT
errorbar(EXPDATA.time{10}(1:4,1), EXPDATA.mean{10}(1:4,1), EXPDATA.SD{10}(1:4,1),'Color', colornormo_lowHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_lowHC, 'MarkerSize',8)
hold on

xticks([0 8 24 48 288 296 312 336 360])
ylim([0 3500])

xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mIU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

subplot(2,2,4)

% Insulin, GTT d13-15

% Hyperglycemia, high HC
errorbar(EXPDATA.time{2}(5:end,1), EXPDATA.mean{2}(5:end,1), EXPDATA.SD{2}(5:end,1),'Color', colorhyper_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_highHC, 'MarkerSize',8)
hold on

% Normoglycemia, high HC, with GTT
errorbar(EXPDATA.time{6}(5:end,1), EXPDATA.mean{6}(5:end,1), EXPDATA.SD{6}(5:end,1),'Color', colornormo_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_highHC, 'MarkerSize',8)
hold on

%Low HC

% Hyperglycemia, low HC
errorbar(EXPDATA.time{8}(5:end,1), EXPDATA.mean{8}(5:end,1), EXPDATA.SD{8}(5:end,1),'Color', colorhyper_lowHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_lowHC, 'MarkerSize',8)
hold on

%  Normoglycemia, low HC, with GTT
errorbar(EXPDATA.time{10}(5:end,1), EXPDATA.mean{10}(5:end,1), EXPDATA.SD{10}(5:end,1),'Color', colornormo_lowHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_lowHC, 'MarkerSize',8)
hold on

% Normoglycemia, low HC, w/o GTT
% errorbar(EXPDATA.time{12}(5:end,1), EXPDATA.mean{12}(5:end,1), EXPDATA.SD{12}(5:end,1),'Color', colornormo_lowHC,...
%     'LineStyle', '--', 'Marker',marker,'MarkerFaceColor',colornormo_lowHC, 'MarkerSize',8)
% hold on

ylim([0 3500])

xticks([0 8 24 48 288 296 312 336 360])

xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mIU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off


% Plot calibration data (low HC)

% Plot glucose data
figure()

title('Calibration data')

% Glucose, GTT d1-3

subplot(2,2,1)

% Hyperglycemia, high HC
errorbar(EXPDATA.time{7}(1:4,1), EXPDATA.mean{7}(1:4,1), EXPDATA.SD{7}(1:4,1),'Color', colorhyper_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_highHC, 'MarkerSize',8)
hold on

xlim([0 48])

subplot(2,2,2)

% Hyperglycemia, high HC
errorbar(EXPDATA.time{7}(5:end,1), EXPDATA.mean{7}(5:end,1), EXPDATA.SD{7}(5:end,1),'Color', colorhyper_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_highHC, 'MarkerSize',8)
hold on

% Normoglycemia, high HC with GTT
errorbar(EXPDATA.time{9}(1:end,1), EXPDATA.mean{9}(1:end,1), EXPDATA.SD{9}(1:end,1),'Color', colornormo_highHC,...
     'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_highHC, 'MarkerSize',8)
hold on

xlim([288 360])  

% Plot insulin data
% Hyperglycemia, high HC

%Insulin, GTT d1-3

subplot(2,2,3)
errorbar(EXPDATA.time{8}(1:4,1), EXPDATA.mean{8}(1:4,1), EXPDATA.SD{8}(1:4,1),'Color', colorhyper_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_highHC, 'MarkerSize',8)
hold on
%errorbar(EXPDATA.time{10}(:,1), EXPDATA.mean{10}(:,1), EXPDATA.SD{10}(:,1),'Color', colornormo_highHC,...
%    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_highHC, 'MarkerSize',8)

xlim([0 48])

% Insulin, GTT d13-15

subplot(2,2,4)

% Hyperglycemia, low HC
errorbar(EXPDATA.time{8}(5:end,1), EXPDATA.mean{8}(5:end,1), EXPDATA.SD{2}(5:end,1),'Color', colorhyper_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colorhyper_highHC, 'MarkerSize',8)
hold on

% Normoglycemia, low HC, with GTT
errorbar(EXPDATA.time{10}(:,1), EXPDATA.mean{10}(:,1), EXPDATA.SD{10}(:,1),'Color', colornormo_highHC,...
    'LineStyle', '-', 'Marker',marker,'MarkerFaceColor',colornormo_highHC, 'MarkerSize',8)
hold on

xlim([288 360])

    