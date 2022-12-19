function plotFunction(param, model,simTime,C)

loadExperimentalData

COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

nTime = length(simTime);
modelname{1} = char(model);

% Simulate
for j = 1 : size(param,1) % For each row of the param set.
    
    try
        for i = 1 : 1
            simData{i} = feval(char(modelname{i}), simTime, C, param(j,:), COSTOPTIONS);        
        end
    catch error
        disp(['Simulation crashed, @ simulation for data: ' num2str(i) ' ... ->' ]);
        error = inf;
        disp(param);
        disp(error)
        return;
    end
    
    plotData.Gmeasured(j,1:nTime) = simData{1}.variablevalues(:,1);
    plotData.Imeasured(j,1:nTime) = simData{1}.variablevalues(:,2);
    plotData.Gliver(j,1:nTime) = simData{1}.variablevalues(:,3);
    plotData.Gislets(j,1:nTime) = simData{1}.variablevalues(:,3);
    plotData.Iliver(j,1:nTime) = simData{1}.variablevalues(:,5);
    plotData.Iislets(j,1:nTime) = simData{1}.variablevalues(:,6);
    
    
end

%% Simulations fitted to experimental measurements
figure('Name', 'Simulations fitted to experimental measurements')

% Glucose in liver compartment
subplot(2,3,1)
hold on
%plot(simTime, plotData.Gliver(1,:), 'r-')
if size(param,1) > 1
    for i = 2 : size(param,1)
        plot(simTime, plotData.Gliver(i,:), 'Color',[0.4940, 0.1840, 0.5560],'LineWidth', 1.25)
    end
end
hold off
title('Liver compartment')
xlabel('Time [h]');
ylabel('Glucose concentration (mM)')
xlim([0 simTime(end)])
set(gca,'FontSize',13)


subplot(2,3,2)
hold on
if size(param,1) > 1
    for i = 2 : size(param,1)
        plot(simTime, plotData.Gliver(i,:), 'Color',[0.4940, 0.1840, 0.5560],'LineWidth', 1.25)
    end
end
hold off
title('Islets compartment')
xlabel('Time [h]');
ylabel('Glucose concentration (mM)')
xlim([0 simTime(end)])
set(gca,'FontSize',13)

% Measured glucose
subplot(2,3,3)
hold on
%plot(simTime, plotData.Gmeasured(1,:), 'r-')
errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1}, 'k.','MarkerSize',14)
if size(param,1) > 1
    for i = 2 : size(param,1)
        plot(simTime, plotData.Gmeasured(i,:), 'Color',[0.4940, 0.1840, 0.5560],'LineWidth', 1.25)
    end
end
hold off
title('Measured')
xlabel('Time [h]');
ylabel('Glucose concentration (mM)')
%legend('Data','Simulation')
xlim([0 simTime(end)])
set(gca,'FontSize',13)


% Insulin liver compartment
subplot (2,3,4)
hold on
%plot(simTime, plotData.Iliver(1,:), 'r-')
if size(param,1) > 1
    for i = 2 : size(param,1)
        plot(simTime, plotData.Iliver(i,:), 'Color',[0.4940, 0.1840, 0.5560],'LineWidth', 1.25)
    end
end
hold off
title('Liver compartment')
xlabel('Time [h]');
ylabel('Insulin concentration (mU/L)')
xlim([0 simTime(end)])
set(gca,'FontSize',13)

% Insulin islets compartment
subplot (2,3,5)
hold on
%plot(simTime, plotData.Iislets(1,:), 'r-')
if size(param,1) > 1
    for i = 2 : size(param,1)
        plot(simTime, plotData.Iislets(i,:), 'Color',[0.4940, 0.1840, 0.5560],'LineWidth', 1.25)
    end
end
hold off
title('Islets compartment')
xlabel('Time [h]');
ylabel('Insulin concentration (mU/L)')
xlim([0 simTime(end)])
set(gca,'FontSize',13)

% Measured insulin
subplot (2,3,6)
hold on
%plot(simTime, plotData.Imeasured(1,:), 'r-')
errorbar(EXPDATA.time{2}, EXPDATA.mean{2}, EXPDATA.SD{2}, 'k.','MarkerSize',14)
if size(param,1) > 1
    for i = 2 : size(param,1)
        plot(simTime, plotData.Imeasured(i,:), 'Color',[0.4940, 0.1840, 0.5560],'LineWidth', 1.25)
    end
end
hold off
title('Measured')
xlabel('Time [h]');
ylabel('Insulin concentration (mU/L)')
legend('Data','Simulation')
xlim([0 simTime(end)])
set(gca,'FontSize',13)

