function plotFunctionUncertainty(param, model, parIndex, param_opt, simTime_11mM, simTime_5p5mM)

loadExperimentalData

COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

%% Plotting settings

color11mM =  [0, 0.4470, 0.7410];         % Blue
color5p5mM = [0.6350, 0.0780, 0.1840];    % Red
color2p8mM = [0.4660, 0.6740, 0.1880];    % Green

marker='o';              % Same marker for all concentrations

% 11 mM glucose

modelname{1} = char(model);
nTime_11mM = length(simTime_11mM);

% Simulate for concentration 11 mM

iC_11mM=[9.519*param(parIndex.i_V_m_hep) 9.519*param(parIndex.i_V_m_islets)...
   131.1*param(parIndex.i_V_m_hep) 131.1*param(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088];

for j = 1 : size(param,1) % For each row of the param set.
    
    try
        for i = 1 : 1
            simData_11mM{i} = feval(char(modelname{i}), simTime_11mM, iC_11mM, param(j,:), COSTOPTIONS);        
        end
    catch error
        disp(['Simulation crashed, @ simulation for data: ' num2str(i) ' ... ->' ]);
        error = inf;
        disp(param);
        disp(error)
        return;
    end
    
    plotData.Gmeasured_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,1);
    plotData.Imeasured_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,2);
    plotData.Gliver_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,3);
    plotData.Gislets_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,4);
    plotData.Iliver_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,5);
    plotData.Iislets_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,6);
    plotData.U_ii_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,27);
    plotData.U_id_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,28);
    plotData.S_i_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,8);
    
    
end

%% Simulate for the optimal parameter values

try
    simData_11mM_opt = feval(char(modelname{i}), simTime_11mM, iC_11mM, param_opt, COSTOPTIONS);        
catch error
    disp(['Simulation crashed, @ simulation for data: ' num2str(i) ' ... ->' ]);
    error = inf;
    disp(param);
    disp(error)
    return;
end

plotData.Gmeasured_11mM_opt = simData_11mM_opt.variablevalues(:,1);
plotData.Imeasured_11mM_opt = simData_11mM_opt.variablevalues(:,2);
plotData.Gliver_11mM_opt = simData_11mM_opt.variablevalues(:,3);
plotData.Gislets_11mM_opt = simData_11mM_opt.variablevalues(:,4);
plotData.Iliver_11mM_opt = simData_11mM_opt.variablevalues(:,5);
plotData.Iislets_11mM_opt = simData_11mM_opt.variablevalues(:,6);
plotData.U_ii_11mM_opt = simData_11mM_opt.variablevalues(:,27);
plotData.U_id_11mM_opt = simData_11mM_opt.variablevalues(:,28);
plotData.S_i_11mM_opt = simData_11mM_opt.variablevalues(:,8);


%% Simulations fitted to experimental measurements

% 11 mM

figure()

subplot(2,3,1)

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gliver_11mM(j,1:nTime_11mM),'Color',color11mM)
    hold on
end

hold off
title('Liver compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xlim([0 simTime_11mM(end)])
ylim([0 12])
set(gca,'TickDir','out','FontSize',20);
box off

% Glucose in the islets compartment
subplot(2,3,2)
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gislets_11mM(j,1:nTime_11mM),'Color',color11mM)
    hold on
end

hold off
title('Islets compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xlim([0 simTime_11mM(end)])
ylim([0 12])
set(gca,'TickDir','out','FontSize',20);
box off

% Measured glucose
subplot(2,3,3)
hold on

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gmeasured_11mM(j,1:nTime_11mM),'Color',color11mM)
    hold on
end

hold on
errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on
plot(simTime_11mM,plotData.Gmeasured_11mM_opt,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime_11mM(end)])
ylim([0 12])
title('Measured')
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

% Insulin in the liver compartment
subplot(2,3,4)
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Iliver_11mM(j,1:nTime_11mM),'Color',color11mM)
    hold on
end

hold off
title('Liver compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Insulin concentration (mU/L)','FontSize',20)
xlim([0 simTime_11mM(end)])
set(gca,'TickDir','out','FontSize',20);
box off

% Insulin in the islets compartment
subplot(2,3,5)
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Iislets_11mM(j,1:nTime_11mM),'Color',color11mM)
    hold on
end

hold off
title('Islets compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Insulin concentration (mU/L)','FontSize',20)
xlim([0 simTime_11mM(end)])
set(gca,'TickDir','out','FontSize',20);
box off

% Measured insulin
subplot(2,3,6)
hold on

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_11mM(j,1:nTime_11mM),'Color',color11mM)
    hold on
end

hold on
errorbar(EXPDATA.time{2}, EXPDATA.mean{2}, EXPDATA.SD{2},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on
plot(simTime_11mM,plotData.Imeasured_11mM_opt,'Color', color11mM,'Linewidth',2);
hold off

xlim([0 simTime_11mM(end)])
title('Measured')
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off


% Insulin resistance

figure()
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.S_i_11mM(j,1:nTime_11mM),'Color',color11mM)
    hold on
end

hold off
title('Islets compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Insulin resistance','FontSize',20)
xlim([0 simTime_11mM(end)])
set(gca,'TickDir','out','FontSize',20);
box off
%% 5.5 mM

nTime_5p5mM = length(simTime_5p5mM);
param(:,parIndex.i_G0)=5.5;
%param(:,parIndex.i_I_GTT)=0;

param_opt(:,parIndex.i_G0)=5.5;
%param_opt(:,parIndex.i_I_GTT)=0;

% Simulate for concentration 5.5 mM

iC_5p5mM=[5.5*param(parIndex.i_V_m_hep) 5.5*param(parIndex.i_V_m_islets)...
   0*param(parIndex.i_V_m_hep) 0*param(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088];%% Calculate maximal and minimal values of the predictions

for j = 1 : size(param,1) % For each row of the param set.
    
    try
        for i = 1 : 1
            simData_5p5mM{i} = feval(char(modelname{i}), simTime_5p5mM, iC_5p5mM, param(j,:), COSTOPTIONS);        
        end
    catch error
        disp(['Simulation crashed, @ simulation for data: ' num2str(i) ' ... ->' ]);
        error = inf;
        disp(param);
        disp(error)
        return;
    end
    
    plotData.Gmeasured_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,1);
    plotData.Imeasured_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,2);
    plotData.Gliver_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,3);
    plotData.Gislets_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,4);
    plotData.Iliver_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,5);
    plotData.Iislets_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,6);
    plotData.U_ii_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,27);
    plotData.U_id_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,28);
    plotData.S_i_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,8);
end

% Simulate for the optimal parameter values

try
    simData_5p5mM_opt = feval(char(modelname{i}), simTime_5p5mM, iC_5p5mM, param_opt, COSTOPTIONS);        
catch error
    disp(['Simulation crashed, @ simulation for data: ' num2str(i) ' ... ->' ]);
    error = inf;
    disp(param);
    disp(error)
    return;
end

plotData.Gmeasured_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,1);
plotData.Imeasured_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,2);
plotData.Gliver_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,3);
plotData.Gislets_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,4);
plotData.Iliver_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,5);
plotData.Iislets_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,6);
plotData.U_ii_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,27);
plotData.U_id_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,28);

% Simulations fitted to experimental measurements
%figure('Name', 'Simulations fitted to experimental measurements')

figure()

% Glucose in liver compartment
subplot(2,3,1)
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gliver_5p5mM(j,1:nTime_11mM),'Color',color5p5mM)
    hold on
end

hold off
title('Liver compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xlim([0 simTime_5p5mM(end)])
ylim([0 12])
set(gca,'TickDir','out','FontSize',20);
box off

subplot(2,3,2)
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gislets_5p5mM(j,1:nTime_11mM),'Color',color5p5mM)
    hold on
end

hold off
title('Islets compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Glucose concentration (mM)','FontSize',20)
xlim([0 simTime_5p5mM(end)])
ylim([0 12])
set(gca,'TickDir','out','FontSize',20);
box off

% Measured glucose
subplot(2,3,3)
hold on

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM)
    hold on
end

hold on
errorbar(EXPDATA.time{5}, EXPDATA.mean{5}, EXPDATA.SD{5},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on
plot(simTime_5p5mM,plotData.Gmeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);
hold off

xlim([0 simTime_5p5mM(end)])
ylim([0 12])
title('Measured')
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

% Insulin liver compartment

subplot(2,3,4)
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Iliver_5p5mM(j,1:nTime_11mM),'Color',color5p5mM)
    hold on
end

hold off
title('Liver compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Insulin concentration (mU/L)','FontSize',20)
xlim([0 simTime_5p5mM(end)])
set(gca,'TickDir','out','FontSize',20);
box off

% Insulin islets compartment
subplot(2,3,5)
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Iislets_5p5mM(j,1:nTime_11mM),'Color',color5p5mM)
    hold on
end

hold off
title('Islets compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Insulin concentration (mU/L)','FontSize',20)
xlim([0 simTime_5p5mM(end)])
set(gca,'TickDir','out','FontSize',20);
box off

% Measured insulin

subplot(2,3,6)
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM)
    hold on
end

hold on
errorbar(EXPDATA.time{6}, EXPDATA.mean{6}, EXPDATA.SD{6},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on
plot(simTime_5p5mM,plotData.Imeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);
hold off

xlim([0 simTime_5p5mM(end)])
title('Measured')
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

% Insulin resistance

figure()
hold on
for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.S_i_5p5mM(j,1:nTime_11mM),'Color',color5p5mM)
    hold on
end

hold off
title('Islets compartment','FontSize',20)
xlabel('Time (h)');
ylabel('Insulin resistance','FontSize',20)
xlim([0 simTime_11mM(end)])
set(gca,'TickDir','out','FontSize',20);
box off

%% Compare glucose consumption between the different concentrations

figure()

% 11 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gmeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    hold on
end

hold on
errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on
plot(simTime_5p5mM,plotData.Gmeasured_11mM_opt,'Color', color11mM,'Linewidth',2);

hold on

% 5.5 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM,'Linewidth',1.5)
    hold on
end

hold on
errorbar(EXPDATA.time{5}, EXPDATA.mean{5}, EXPDATA.SD{5},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on
plot(simTime_5p5mM,plotData.Gmeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);

xlim([0 287.9])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off


%% Compare insulin consumption between the different concentrations

figure()

% 11 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    hold on
end

hold on
errorbar(EXPDATA.time{2}, EXPDATA.mean{2}, EXPDATA.SD{2},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on
plot(simTime_5p5mM,plotData.Imeasured_11mM_opt,'Color', color11mM,'Linewidth',2);

hold on

% 5.5 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM,'Linewidth',1.5)
    hold on
end

hold on
errorbar(EXPDATA.time{6}, EXPDATA.mean{6}, EXPDATA.SD{6},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on
plot(simTime_5p5mM,plotData.Imeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);

xlim([0 287.9])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

%% Compare glucose consumption between the different concentrations for the experiment with dose A

figure()

% 11 mM

for j = 1 : size(param,1)
    f1=plot(simTime_11mM,plotData.Gmeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    f1.Color(4) = 0.25;
    hold on
end

hold on

plot(simTime_11mM,plotData.Gmeasured_11mM_opt,'Color', color11mM,'Linewidth',2);

hold on

% 5.5 mM

for j = 1 : size(param,1)
    f2=plot(simTime_11mM,plotData.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM,'Linewidth',1.5)
    f2.Color(4) = 0.25;
    hold on
end

plot(simTime_5p5mM,plotData.Gmeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);

xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

% Plot experiment data

hold on
errorbar(EXPDATA.time{7}, EXPDATA.mean{7}, EXPDATA.SD{7},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on
errorbar(EXPDATA.time{8}, EXPDATA.mean{8}, EXPDATA.SD{8},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)

%% Compare glucose consumption between the different concentrations for the experiment with dose B

figure()

% 11 mM

for j = 1 : size(param,1)
    f1=plot(simTime_11mM,plotData.Gmeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    f1.Color(4) = 0.25;
    hold on
end

hold on

plot(simTime_11mM,plotData.Gmeasured_11mM_opt,'Color', color11mM,'Linewidth',2);

hold on

% 5.5 mM

for j = 1 : size(param,1)
    f2=plot(simTime_11mM,plotData.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM,'Linewidth',1.5)
    f2.Color(4) = 0.25;
    hold on
end

plot(simTime_5p5mM,plotData.Gmeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);

xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

% Plot experiment data

hold on
errorbar(EXPDATA.time{9}, EXPDATA.mean{9}, EXPDATA.SD{9},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on
errorbar(EXPDATA.time{10}, EXPDATA.mean{10}, EXPDATA.SD{10},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)%% Compare insulin consumption between the different concentrations for the experiment


%% Compare insulin consumption between the different concentrations for the experiment

figure()

% 11 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    hold on
end

hold on

plot(simTime_11mM,plotData.Imeasured_11mM_opt,'Color', color11mM,'Linewidth',2,'Linewidth',2);

hold on

% 5.5 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_5p5mM(j,1:nTime_5p5mM),'Color',color5p5mM,'Linewidth',1.5)
    hold on
end

hold on
plot(simTime_5p5mM,plotData.Imeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);

hold on

xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

% Plot experiment data

hold on
errorbar(EXPDATA.time{13}, EXPDATA.mean{13}, EXPDATA.SD{13},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on
errorbar(EXPDATA.time{14}, EXPDATA.mean{14}, EXPDATA.SD{14},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)

figure()

% 11 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    hold on
end

hold on

plot(simTime_11mM,plotData.Imeasured_11mM_opt,'Color', color11mM,'Linewidth',2);

hold on

% 5.5 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_5p5mM(j,1:nTime_5p5mM),'Color',color5p5mM,'Linewidth',1.5)
    hold on
end

hold on
plot(simTime_5p5mM,plotData.Imeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);

hold on

xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

%% Compare glucose consumption between the different concentrations for the experiment

figure()

% 11 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gmeasured_11mM(j,1:nTime_11mM)-plotData.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color2p8mM,'LineWidth',1.5)
    hold on
end
hold on

% plot(simTime_11mM,plotData.Gmeasured_11mM_opt,'Color', color11mM,'Linewidth',2);

hold on

% 5.5 mM

xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

%% Compare insulin consumption between the different concentrations for the experiment

figure()

% 11 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    hold on
end

hold on

plot(simTime_11mM,plotData.Imeasured_11mM_opt,'Color', color11mM,'Linewidth',2,'Linewidth',2);

hold on

% 5.5 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Imeasured_5p5mM(j,1:nTime_5p5mM),'Color',color5p5mM,'Linewidth',1.5)
    hold on
end

hold on
plot(simTime_5p5mM,plotData.Imeasured_5p5mM_opt,'Color', color5p5mM,'Linewidth',2);

hold on

xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

