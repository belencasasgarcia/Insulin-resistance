function plotData = plotFunctionExperimentH1(param, model, parIndex, param_opt, simTime_11mM, simTime_5p5mM, dose, limit_iiuptake,G_GTT_hyper,G_GTT_normo,N_samples_G0,...
    SEM_G_hyper,SEM_G_normo)

loadExperimentalData

COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;


%% Plotting settings

color11mM =  [0.6350, 0.0780, 0.1840];    % Red
color5p5mM = [0.4660, 0.6740, 0.1880];    % Green

marker='o';              % Same marker for all concentrations

% 11 mM glucose

modelname{1} = char(model);
nTime_11mM = length(simTime_11mM);

param(:,parIndex.i_G_GTT)=G_GTT_hyper;
param(:,parIndex.i_delta_G_7)=param(:,parIndex.i_delta_G_7_11mM);
param(:,parIndex.i_delta_I_7)=param(:,parIndex.i_delta_I_7_11mM);
param(:,parIndex.i_I_GTT_13)=dose;

param_IC=repmat(param,N_samples_G0,1); %Sample N_samples values of glucose concentrations
% between (mean-SEM,mean+SEM)

param_IC0=[linspace(-SEM_G_hyper,0,N_samples_G0/2) linspace(0,SEM_G_hyper,N_samples_G0/2)];
param_IC0=repmat(param_IC0,size(param,1),1);
param_IC0=reshape(param_IC0,[],1);


param_IC(:,parIndex.i_G_GTT)=repmat(G_GTT_hyper,size(param,1)*N_samples_G0,1)+param_IC0;

% Simulate for concentration 11 mM

for j = 1 : size(param_IC,1) % For each row of the param set.
    
    iC_11mM=[(11+param_IC(j,parIndex.i_delta_G_1_11mM))*param_IC(j,parIndex.i_V_m_hep) ...
            (11+param_IC(j,parIndex.i_delta_G_1_11mM))*param_IC(j,parIndex.i_V_m_islets)...
            (0+param_IC(j,parIndex.i_delta_I_1_11mM))*param_IC(j,parIndex.i_V_m_hep) ...
            (0+param_IC(j,parIndex.i_delta_I_1_11mM))*param_IC(j,parIndex.i_V_m_islets) 0 0 ...
            5.5 0.0000000088 0 0];
    try
        for i = 1 : 1
            simData_11mM{i} = feval(char(modelname{i}), simTime_11mM, iC_11mM, param_IC(j,:), COSTOPTIONS);        
        end
    catch error
        disp(['Simulation crashed, @ simulation for data: ' num2str(i) ' ... ->' ]);
        error = inf;
        disp(param);
        disp(error)
        return;
    end
    
    plotData.Guptake_id_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,31);
    plotData.Guptake_ii_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,32);
    plotData.Guptake_ratio_11mM(j,1:nTime_11mM)=plotData.Guptake_id_11mM(j,1:nTime_11mM)./...
        (plotData.Guptake_ii_11mM(j,1:nTime_11mM)+plotData.Guptake_id_11mM(j,1:nTime_11mM))*100;
end

% Find parameter values that have insulin dependent uptake> 50% of total
% uptake during at the end of the first GTT

ratio_48h=plotData.Guptake_ratio_11mM(:,4800); % Uptake ratio at the end of the first GTT
indexes=ratio_48h>=limit_iiuptake;
param_IC=param_IC(indexes,:);

for j = 1 : size(param,1) % For each row of the param set.
    
    iC_11mM=[(11+param(j,parIndex.i_delta_G_1_11mM))*param(j,parIndex.i_V_m_hep) ...
            (11+param(j,parIndex.i_delta_G_1_11mM))*param(j,parIndex.i_V_m_islets)...
            (0+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_hep) ...
            (0+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_islets) 0 0 ...
            5.5 0.0000000088 0 0];
    try
        for i = 1 : 1
            simData_11mM{i} = feval(char(modelname{i}), simTime_11mM, iC_11mM, param_IC(j,:), COSTOPTIONS);        
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

%% 5.5 mM

nTime_5p5mM = length(simTime_5p5mM);

param_IC(:,parIndex.i_G_GTT)=G_GTT_normo;

param_IC(:,parIndex.i_G0)=5.5;
param_IC(:,parIndex.i_delta_G_7)=-0.315;

% Simulate for concentration 5.5 mM

iC_5p5mM=[5.5*param_IC(parIndex.i_V_m_hep) 5.5*param_IC(parIndex.i_V_m_islets)...
   0*param_IC(parIndex.i_V_m_hep) 0*param_IC(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];%% Calculate maximal and minimal values of the predictions

for j = 1 : size(param,1) % For each row of the param set.
    
    try
        for i = 1 : 1
            simData_5p5mM{i} = feval(char(modelname{i}), simTime_5p5mM, iC_5p5mM, param_IC(j,:), COSTOPTIONS);        
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

% Plot prediction difference between normo and hyperglycemia (for dose
% calculation)

figure()

% 11 mM

for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gmeasured_11mM(j,1:nTime_11mM)-plotData.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    hold on
end

xlim([288 simTime_11mM(end)])
ylim([0 13.75])
yticks([0 2.75 5.5 11 13.75])
xlabel('Time (h)','FontSize',20);
set(gca,'TickDir','out','FontSize',20);
box off

hold on


