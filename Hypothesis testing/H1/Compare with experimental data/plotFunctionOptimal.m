function plotFunctionOptimal(param, model)

loadExperimentalData


COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

simTime = [0:0.1:48];
nTime = length(simTime);
modelname{1} = char(model);

iCOrig = [11*param(1) 11*param(3) 0 0];

% Simulate
for j = 1 : size(param,1) % For each row of the param set.
    
    try
        for i = 1 : 1
            simData{i} = feval(char(modelname{i}), simTime, iCOrig, param(j,:), COSTOPTIONS);        
        end
    catch error
        disp(['Simulation crashed, @ simulation for data: ' num2str(i) ' ... ->' ]);
        error = inf;
        disp(param);
        disp(error)
        return;
    end
    
    
    plotData.Gliver(j,1:nTime) = simData{1}.variablevalues(:,1);
    plotData.Gislets(j,1:nTime) = simData{1}.variablevalues(:,3);
    plotData.Iliver(j,1:nTime) = simData{1}.variablevalues(:,4);
    plotData.Iislets(j,1:nTime) = simData{1}.variablevalues(:,2);
    plotData.Gmeasured(j,1:nTime) = simData{1}.variablevalues(:,5);
    plotData.Imeasured(j,1:nTime) = simData{1}.variablevalues(:,6);
    
end

%% Simulations fitted to experimental measurements
figure('Name', 'Simulations fitted to experimental measurements')

% Glucose in liver compartment
subplot(2,3,1)
hold on
plot(simTime, plotData.Gliver(1,:), 'b-')
hold off
title('Liver compartment')
xlabel('Time [h]');
ylabel('Glucose concentration (mM)')
legend('Simulated')
xlim([0 simTime(end)])


subplot(2,3,2)
hold on
plot(simTime, plotData.Gislets(1,:), 'b-')
hold off
title('Islets compartment')
xlabel('Time [h]');
ylabel('Glucose concentration (mM)')
legend('Simulated')
xlim([0 simTime(end)])

% Measured glucose
subplot(2,3,3)
hold on
plot(simTime, plotData.Gmeasured(1,:), 'b-')
errorbar(EXPDATA.time{1}, EXPDATA.mean{1}, EXPDATA.SD{1}, 'k^--')
hold off
title('Measured')
xlabel('Time [h]');
ylabel('Glucose concentration (mM)')
legend('Measured','Simulated')
xlim([0 simTime(end)])


% Insulin liver compartment
subplot (2,3,4)
hold on
plot(simTime, plotData.Iliver(1,:), 'b-')
hold off
title('Liver compartment')
xlabel('Time [h]');
ylabel('Insulin concentration (mU/L)')
legend('Simulated')
xlim([0 simTime(end)])

% Insulin islets compartment
subplot (2,3,5)
hold on
plot(simTime, plotData.Iislets(1,:), 'b-')
hold off
title('Islets compartment')
xlabel('Time [h]');
ylabel('Insulin concentration (mU/L)')
legend('Simulated')
xlim([0 simTime(end)])

% Measured insulin
subplot (2,3,6)
hold on
plot(simTime, plotData.Imeasured(1,:), 'b-')
errorbar(EXPDATA.time{2}, EXPDATA.mean{2}, EXPDATA.SD{2}, 'k^--')
hold off
title('Measured')
xlabel('Time [h]');
ylabel('Insulin concentration (mU/L)')
legend('Measured','Simulated')
xlim([0 simTime(end)])

% Glucose, only hepatocytes present


