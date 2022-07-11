function plotFunctionUncertainty(simTime_11mM, simTime_5p5mM, plotData_H1, plotData_H2)

%% Plotting settings

color11mM =  [0, 0.4470, 0.7410];         % Blue
color5p5mM = [0.6350, 0.0780, 0.1840];    % Red
color_difH1 = [0.4660, 0.6740, 0.1880];   % Green
color_difH2 = [0.8500, 0.3250, 0.0980];   % Orange

marker='o';              % Same marker for all concentrations

nTime_11mM = length(simTime_11mM);
nTime_5p5mM = length(simTime_5p5mM);

%% Compare glucose consumption between the different concentrations for the experiment

figure()

%H1

for j = 1 : 60
    plot(simTime_11mM,plotData_H1.Gmeasured_11mM(j,1:nTime_11mM)-plotData_H1.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color_difH1,'LineWidth',1.5)
    hold on
end

hold on

%H2

for j = 1 : 60
    plot(simTime_11mM,plotData_H2.Gmeasured_11mM(j,1:nTime_11mM)-plotData_H2.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color_difH2,'LineWidth',1.5)
    hold on
end


xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

