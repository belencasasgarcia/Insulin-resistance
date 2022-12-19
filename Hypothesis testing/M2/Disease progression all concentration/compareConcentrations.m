function [] = compareConcentrations(time1,data1_max,data1_min,EXPDATA_index1,...
    time2,data2_max,data2_min,EXPDATA_index2,time3,data3_max,data3_min,EXPDATA_index3)

color11mM =  [0, 0.4470, 0.7410];         % Blue
color5p5mM = [0.6350, 0.0780, 0.1840];    % Red
color2p8mM = [0.4660, 0.6740, 0.1880];    % Green

marker='o';              % Same marker for all concentrations

%% Compare insulin consumption between the different concentrations

figure()

% 11 mM

area1=[time1,fliplr(time1)];
area2=[data1_max,fliplr(data1_min)];
f=fill(area1,area2,color11mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{EXPDATA_index1}, EXPDATA.mean{EXPDATA_index1}, ...
    EXPDATA.SD{EXPDATA_index1},'Color', color11mM, 'LineStyle', 'n', ...
    'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on
plot(simTime_5p5mM,plotData.Imeasured_11mM_opt,'Color', color11mM,'Linewidth',2);

hold on

% 5.5 mM

area1=[simTime_5p5mM,fliplr(simTime_5p5mM)];
area2=[plotData.Imeasured_5p5mM_max,fliplr(plotData.Gmeasured_5p5mM_min)];
f=fill(area1,area2,color5p5mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{5}, EXPDATA.mean{5}, EXPDATA.SD{5},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on
plot(simTime_5p5mM,plotData.Gmeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);

hold on

area1=[simTime_2p8mM,fliplr(simTime_2p8mM)];
area2=[plotData.Gmeasured_2p8mM_max,fliplr(plotData.Gmeasured_2p8mM_min)];
f=fill(area1,area2,color2p8mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{7}, EXPDATA.mean{7}, EXPDATA.SD{7},'Color', color2p8mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color2p8mM, 'MarkerSize',8)
hold on
plot(simTime_2p8mM,plotData.Gmeasured_2p8mM_opt,'Color', color2p8mM,'Linewidth',2);
hold off

xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off
end

