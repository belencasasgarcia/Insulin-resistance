function plotFunctionUncertainty_lowHC(param, model, parIndex, param_opt, simTime_11mM, simTime_5p5mM)

loadExperimentalData_HC
time_end = 336;

COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

%% Plotting settings

color11mM =  [0.6350, 0.0780, 0.1840];    % Blue
color5p5mM = [0.4660, 0.6740, 0.1880];    % Green

marker='o';              % Same marker for all concentrations

% 11 mM glucose

modelname{1} = char(model);
nTime_11mM = length(simTime_11mM);

offset_11mM = 0
offset_5p5mM = 0

% Simulate for concentration 11 mM

for j = 1 : size(param,1) % For each row of the param set.

    % Change insulin secretion capacity and insulin sensitivity to match
    % optimal values from the optimization corresponding to the
    % optimization from GTT d1-3, low HCT. 

    % The parameters correspond to one parameter set from the ...
    % 'opt_par_GTTd1to3_spreaded.mat' dataset

    param(j,parIndex.i_Sigma)=1.5e+07;
    param(j,parIndex.i_S_i)=5.67e-04;
    
    param(j,parIndex.i_delta_G_1_11mM)=0;
    param(j,parIndex.i_delta_G_1_11mM)=0;
    param(j,parIndex.i_delta_G_13_11mM)=offset_11mM;
    param(j,parIndex.i_delta_I_1_11mM)=0;

    iC_11mM=[(11+param(j,parIndex.i_delta_G_1_11mM))*...
        param(j,parIndex.i_V_m_hep) (11+param(j,parIndex.i_delta_G_1_11mM))...
        *param(j,parIndex.i_V_m_islets)...
   (11+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_hep) ...
   (11+param(j,parIndex.i_delta_I_1_11mM))*param(j,parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];

   param(j,parIndex.i_delta_G_13) = param(j,parIndex.i_delta_G_13_11mM);
   param(j,parIndex.i_delta_I_13) = param(j,parIndex.i_delta_I_13_11mM);

    try
        for i = 1 : 1
            simData_11mM{i} = feval(char(modelname{i}), simTime_11mM, ...
                iC_11mM, param(j,:), COSTOPTIONS);        
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
    plotData.ratio_id_ii_11mM(j,1:nTime_11mM)=plotData.U_id_11mM(j,1:nTime_11mM)...
        ./plotData.U_ii_11mM(j,1:nTime_11mM);
    

end

%% Calculate maximal and minimal values of the predictions

plotData.Gmeasured_11mM_max(1:nTime_11mM)=max(plotData.Gmeasured_11mM,[],1);
plotData.Gmeasured_11mM_min(1:nTime_11mM)=min(plotData.Gmeasured_11mM,[],1);
plotData.Imeasured_11mM_max(1:nTime_11mM)=max(plotData.Imeasured_11mM,[],1);
plotData.Imeasured_11mM_min(1:nTime_11mM)=min(plotData.Imeasured_11mM,[],1);
plotData.Gliver_11mM_max(1:nTime_11mM)=max(plotData.Gliver_11mM,[],1);
plotData.Gliver_11mM_min(1:nTime_11mM)=min(plotData.Gliver_11mM,[],1);
plotData.Gislets_11mM_max(1:nTime_11mM)=max(plotData.Gislets_11mM,[],1);
plotData.Gislets_11mM_min(1:nTime_11mM)=min(plotData.Gislets_11mM,[],1);
plotData.Iliver_11mM_max(1:nTime_11mM)=max(plotData.Iliver_11mM,[],1);
plotData.Iliver_11mM_min(1:nTime_11mM)=min(plotData.Iliver_11mM,[],1);
plotData.Iislets_11mM_max(1:nTime_11mM)=max(plotData.Iislets_11mM,[],1);
plotData.Iislets_11mM_min(1:nTime_11mM)=min(plotData.Iislets_11mM,[],1);
plotData.U_ii_11mM_min(1:nTime_11mM)=min(plotData.U_ii_11mM,[],1);
plotData.U_ii_11mM_max(1:nTime_11mM)=max(plotData.U_ii_11mM,[],1);
plotData.U_id_11mM_min(1:nTime_11mM)=min(plotData.U_id_11mM,[],1);
plotData.U_id_11mM_max(1:nTime_11mM)=max(plotData.U_id_11mM,[],1);
plotData.S_i_11mM_max(j,1:nTime_11mM) = max(plotData.S_i_11mM,[],1);
plotData.S_i_11mM_min(j,1:nTime_11mM) = min(plotData.S_i_11mM,[],1);
plotData.G_slow_11mM_max(j,1:nTime_11mM) = max(plotData.G_slow_11mM,[],1);
plotData.G_slow_11mM_min(j,1:nTime_11mM) = min(plotData.G_slow_11mM,[],1);
plotData.beta_vol_11mM_max(j,1:nTime_11mM) = max(plotData.beta_vol_11mM,[],1);
plotData.beta_vol_11mM_min(j,1:nTime_11mM) = min(plotData.beta_vol_11mM,[],1);
plotData.S_cap_11mM_max(j,1:nTime_11mM) = max(plotData.S_cap_11mM,[],1);
plotData.S_cap_11mM_min(j,1:nTime_11mM) = min(plotData.S_cap_11mM,[],1);
plotData.S_cap_11mM_max(j,1:nTime_11mM) = max(plotData.S_cap_11mM,[],1);
plotData.ratio_id_ii_11mM_min(j,1:nTime_11mM) = min(plotData.ratio_id_ii_11mM,[],1);
plotData.ratio_id_ii_11mM_max(j,1:nTime_11mM) = max(plotData.ratio_id_ii_11mM,[],1);
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
plotData.G_slow_opt = simData_11mM_opt.variablevalues(:,9);
plotData.beta_vol_opt = simData_11mM_opt.variablevalues(:,12);
plotData.S_cap_opt = simData_11mM_opt.variablevalues(:,14);

%% 5.5 mM

nTime_5p5mM = length(simTime_5p5mM);
param(:,parIndex.i_G0)=5.5;
param_opt(:,parIndex.i_G0)=5.5;

% Simulate for concentration 5.5 mM

for j = 1 : size(param,1) % For each row of the param set.

    param(j,parIndex.i_delta_G_1_5p5mM)=0;
    param(j,parIndex.i_delta_I_1_5p5mM)=0;

    param(j,parIndex.i_delta_G_13_5p5mM)=offset_5p5mM;
    param(j,parIndex.i_delta_I_13_5p5mM)=0;


    iC_5p5mM=[(5.5+param(j,parIndex.i_delta_G_1_5p5mM))*param(j,parIndex.i_V_m_hep) (5.5+param(j,parIndex.i_delta_G_1_5p5mM))*param(j,parIndex.i_V_m_islets)...
   (0+param(j,parIndex.i_delta_I_1_5p5mM))*param(j,parIndex.i_V_m_hep) (0+param(j,parIndex.i_delta_I_1_5p5mM))*param(j,parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];

   param(:,parIndex.i_delta_G_13)=param(:,parIndex.i_delta_G_13_5p5mM);
   param(:,parIndex.i_delta_I_13)=param(:,parIndex.i_delta_I_13_5p5mM);
    
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
    plotData.int_G_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,7);
    plotData.S_i_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,8);
    plotData.G_slow_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,9);
    plotData.beta_vol_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,12);
    plotData.S_cap_5p5mM(j,1:nTime_11mM) = simData_5p5mM{1}.variablevalues(:,14);
end

plotData.Gmeasured_5p5mM_max(1:nTime_5p5mM)=max(plotData.Gmeasured_5p5mM,[],1);
plotData.Gmeasured_5p5mM_min(1:nTime_5p5mM)=min(plotData.Gmeasured_5p5mM,[],1);
plotData.Imeasured_5p5mM_max(1:nTime_5p5mM)=max(plotData.Imeasured_5p5mM,[],1);
plotData.Imeasured_5p5mM_min(1:nTime_5p5mM)=min(plotData.Imeasured_5p5mM,[],1);
plotData.Gliver_5p5mM_max(1:nTime_5p5mM)=max(plotData.Gliver_5p5mM,[],1);
plotData.Gliver_5p5mM_min(1:nTime_5p5mM)=min(plotData.Gliver_5p5mM,[],1);
plotData.Gislets_5p5mM_max(1:nTime_5p5mM)=max(plotData.Gislets_5p5mM,[],1);
plotData.Gislets_5p5mM_min(1:nTime_5p5mM)=min(plotData.Gislets_5p5mM,[],1);
plotData.Iliver_5p5mM_max(1:nTime_5p5mM)=max(plotData.Iliver_5p5mM,[],1);
plotData.Iliver_5p5mM_min(1:nTime_5p5mM)=min(plotData.Iliver_5p5mM,[],1);
plotData.Iislets_5p5mM_max(1:nTime_5p5mM)=max(plotData.Iislets_5p5mM,[],1);
plotData.Iislets_5p5mM_min(1:nTime_5p5mM)=min(plotData.Iislets_5p5mM,[],1);
plotData.U_ii_5p5mM_min(1:nTime_11mM)=min(plotData.U_ii_5p5mM,[],1);
plotData.U_ii_5p5mM_max(1:nTime_11mM)=max(plotData.U_ii_5p5mM,[],1);
plotData.U_id_5p5mM_min(1:nTime_11mM)=min(plotData.U_id_5p5mM,[],1);
plotData.U_id_5p5mM_max(1:nTime_11mM)=max(plotData.U_id_5p5mM,[],1);
plotData.int_G_5p5mM_max(j,1:nTime_11mM) = max(plotData.int_G_5p5mM,[],1);
plotData.int_G_5p5mM_min(j,1:nTime_11mM) = min(plotData.int_G_5p5mM,[],1);
plotData.S_i_5p5mM_max(j,1:nTime_11mM) = max(plotData.S_i_5p5mM,[],1);
plotData.S_i_5p5mM_min(j,1:nTime_11mM) = min(plotData.S_i_5p5mM,[],1);
plotData.G_slow_5p5mM_max(j,1:nTime_11mM) = max(plotData.G_slow_5p5mM,[],1);
plotData.G_slow_5p5mM_min(j,1:nTime_11mM) = min(plotData.G_slow_5p5mM,[],1);
plotData.beta_vol_5p5mM_max(j,1:nTime_11mM) = max(plotData.beta_vol_5p5mM,[],1);
plotData.beta_vol_5p5mM_min(j,1:nTime_11mM) = min(plotData.beta_vol_5p5mM,[],1);
plotData.S_cap_5p5mM_max(j,1:nTime_11mM) = max(plotData.S_cap_5p5mM,[],1);
plotData.S_cap_5p5mM_min(j,1:nTime_11mM) = min(plotData.S_cap_5p5mM,[],1);


% Simulate for the optimal parameter values

param_opt(:,parIndex.i_delta_G_13) = param_opt(:,parIndex.i_delta_G_13_5p5mM);
param_opt(:,parIndex.i_delta_I_13) = param_opt(:,parIndex.i_delta_I_13_5p5mM);

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
plotData.S_i_5p5mM_opt = simData_5p5mM_opt.variablevalues(:,8);
plotData.G_slow_5p5mM_opt = simData_11mM_opt.variablevalues(:,9);
plotData.beta_vol_5p5mM_opt = simData_11mM_opt.variablevalues(:,12);
plotData.S_cap_5p5mM_opt = simData_11mM_opt.variablevalues(:,14);

%% Plot GTTs on days 1 and 13
%% Compare glucose consumption between the different concentrations

figure()
subplot(1,2,1)

% Plot GTT day 1-3

% 11 mM

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[plotData.Gmeasured_11mM_max,fliplr(plotData.Gmeasured_11mM_min)];
f=fill(area1,area2,color11mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{7}, EXPDATA.mean{7}, EXPDATA.SD{7},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on

% 5.5 mM

xlim([0 48])
ylim([0 16])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

% Plot GTT day 13-15

% 11 mM
subplot(1,2,2)

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[plotData.Gmeasured_11mM_max,fliplr(plotData.Gmeasured_11mM_min)];
f=fill(area1,area2,color11mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{7}, EXPDATA.mean{7}, EXPDATA.SD{7},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on

% 5.5 mM

area1=[simTime_5p5mM,fliplr(simTime_5p5mM)];
area2=[plotData.Gmeasured_5p5mM_max,fliplr(plotData.Gmeasured_5p5mM_min)];
f=fill(area1,area2,color5p5mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{9}, EXPDATA.mean{9}, EXPDATA.SD{9},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)

xlim([288 time_end])
ylim([0 16])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose concentration (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

%% Compare insulin consumption between the different concentrations

figure()

subplot(1,2,1)
% Plot GTT day 13-15

% 11 mM

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[plotData.Imeasured_11mM_max,fliplr(plotData.Imeasured_11mM_min)];
f=fill(area1,area2,color11mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{8}, EXPDATA.mean{8}, EXPDATA.SD{8},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on

% 5.5 mM

xlim([0 48])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mIU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

subplot(1,2,2)
% Plot GTT day 13-15
% 11 mM

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[plotData.Imeasured_11mM_max,fliplr(plotData.Imeasured_11mM_min)];
f=fill(area1,area2,color11mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{8}, EXPDATA.mean{8}, EXPDATA.SD{8},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)

hold on

% 5.5 mM

area1=[simTime_5p5mM,fliplr(simTime_5p5mM)];
area2=[plotData.Imeasured_5p5mM_max,fliplr(plotData.Imeasured_5p5mM_min)];
f=fill(area1,area2,color5p5mM);
alpha(f,.2)

hold on
errorbar(EXPDATA.time{10}, EXPDATA.mean{10}, EXPDATA.SD{10},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)

xlim([0 time_end])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off

xlim([288 time_end])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin concentration (mIU/L)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
box off