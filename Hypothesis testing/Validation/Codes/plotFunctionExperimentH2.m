function plotData = plotFunctionExperimentH2 (param, model, parIndex, param_opt, simTime_11mM, simTime_5p5mM,dose,limit_iiuptake,ratio_F_hyper,G_GTT_hyper,G_GTT_normo)

%   plotFunctionExperimentH1 Plots prediction of glucose responses for
%   hypothesis H2 given an insulin spike "dose"

%   Input parameters:
%       param: Acceptable parameter values
%       model: Model to simulate (hypothesis H1)
%       param_Index: Parameter indexes in the model
%       param_opt: Optimal parameter values
%       simTime_11mM: Time axis to simulate hyperglycemic concentrations
%       simTime_5p5mM: Time axis to simulate normo concentrations
%       dose: Insulin dose (mIU/L)
%       limit_iiuptake: Limit in insulin independent glucose uptake for the
%       prediction
%       G_GTT_hyper: Glucose value for GTT under hyperglycemia 
%       G_GTT_normo: Glucose value for GTT under normoglycemia 
%       N_samples_G0: Number of samples to simulate between (-SEM,SEM)
%       SEM_G_hyper: SEM in experimental data for hyperglycemia
%       SEM_G_normo: SEM in experimental data for normoglycemia

%   Output:
%       plotData: Simulated data

% Define options to simulate the model
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

% Simulate for concentration 11 mM

for j = 1 : size(param,1) % For each row of the param set.
    
    iC_11mM=[(11+param(j,parIndex.i_delta_G_1_11mM))*param(j,parIndex.i_V_m_hep) ...
            (11+param(j,parIndex.i_delta_G_1_11mM))*param(j,parIndex.i_V_m_islets)...
            (0+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_hep) ...
            (0+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_islets) 0 0 ...
            5.5 0.0000000088 0 0];
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
    
    plotData.Guptake_id_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,31);
    plotData.Guptake_ii_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,32);
    plotData.Guptake_ratio_11mM(j,1:nTime_11mM)=plotData.Guptake_id_11mM(j,1:nTime_11mM)./...
        (plotData.Guptake_ii_11mM(j,1:nTime_11mM)+plotData.Guptake_id_11mM(j,1:nTime_11mM))*100;
    plotData.F_hyper_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,29);

end

ratio_48h=plotData.Guptake_ratio_11mM(:,4800); % Uptake ratio at the end of the first GTT
indexes=plotData.F_hyper_11mM(:,end)>ratio_F_hyper & ratio_48h>=limit_iiuptake;
param=param(indexes,:);


for j = 1 : size(param,1) % For each row of the param set.
    
    iC_11mM=[(11+param(j,parIndex.i_delta_G_1_11mM))*param(j,parIndex.i_V_m_hep) ...
            (11+param(j,parIndex.i_delta_G_1_11mM))*param(j,parIndex.i_V_m_islets)...
            (0+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_hep) ...
            (0+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_islets) 0 0 ...
            5.5 0.0000000088 0 0];
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
    plotData.G_slow_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,9);
    plotData.beta_vol_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,12);
    plotData.S_cap_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,14);
    plotData.Guptake_id_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,31);
    plotData.Guptake_ii_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,32);
    plotData.Guptake_ratio_11mM(j,1:nTime_11mM)=plotData.Guptake_id_11mM(j,1:nTime_11mM)./...
        (plotData.Guptake_ii_11mM(j,1:nTime_11mM)+plotData.Guptake_id_11mM(j,1:nTime_11mM))*100;
    plotData.F_hyper_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,29);
    plotData.F_additional_11mM(j,1:nTime_11mM) = simData_11mM{1}.variablevalues(:,30);
    

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


%% 5.5 mM

nTime_5p5mM = length(simTime_5p5mM);
param(:,parIndex.i_G_GTT)=G_GTT_normo;
param(:,parIndex.i_G0)=5.5;
param(:,parIndex.i_delta_G_7)=-0.315;

% Simulate for concentration 5.5 mM

iC_5p5mM=[5.5*param(parIndex.i_V_m_hep) 5.5*param(parIndex.i_V_m_islets)...
   0*param(parIndex.i_V_m_hep) 0*param(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];%% Calculate maximal and minimal values of the predictions

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
    plotData.G_slow_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,9);
    plotData.beta_vol_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,12);
    plotData.S_cap_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,14);
    plotData.Guptake_id_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,31);
    plotData.Guptake_ii_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,32);
    plotData.Guptake_ratio_5p5mM(j,1:nTime_11mM)=plotData.Guptake_id_5p5mM(j,1:nTime_11mM)./(plotData.Guptake_ii_11mM(j,1:nTime_11mM)...
        +plotData.Guptake_id_11mM(j,1:nTime_11mM))*100;
    plotData.F_hyper_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,29);
    plotData.F_additional_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,30);
    
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


for j = 1 : size(param,1)
    plot(simTime_11mM,plotData.Gmeasured_11mM(j,1:nTime_11mM)-plotData.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM,'Linewidth',1.5)
    hold on
end

xlim([288 336])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose (mM)','FontSize',20);
set(gca,'TickDir','out','FontSize',20);
box off


