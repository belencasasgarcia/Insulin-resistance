clear
clc
close all
format long
format compact

% Declare global variables

global VALIDATIONDATA
global Hypothesis1
global Hypothesis2
global icOrig0
global pNames
global param
global OPTIONS
global COSTOPTIONS
global EXPDATA
global SIMTIME      % Time axis for simulation
global parIndex

%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load experimental data

loadExperimentalData

k_ins_nM=1/144;

time_end=360;

SIMTIME=[0:0.005:time_end];

limit_iiuptake=45;
ratio_F_hyper=0.6;

%% Parameter indexes

% Load parameter indexes for the model corresponding to hypothesis 1
loadParameterIndexesH1

% Start conditions
iC=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
   183.3*param(parIndex.i_V_m_hep) 183.3*param(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];

% SBPDsimulate options
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through

% Settings for plots
marker='o';  

color11mM = [204, 0, 51]/255;                        % Hyperglycemia, 50 µM cortisone
colorhyper_lowHC = [0, 64, 128]/255;                 % Hyperglycemia, 10 nM cortisone
colornormo_lowHC = [149, 214, 246]/255;              % Normoglycemia, 10 nM cortisone
color5p5mM = [235, 174, 9]/255;                      % Normoglycemia, 50 µM cortisone

colorH1 = [0, 0.4470, 0.7410];                       % Blue
colorH2 = [0.8500, 0.3250, 0.0980];                  % Orange

%%  Dose 1

% Insulin dose value
% Simulations were performed for two insulin doses: dose A (2500 mIU/L) 
% and dose B (3500 mIU/L). In the publication, only the second dose is included

dose = [2500 3500]; % mIU/L

% Start glucose values for the GTT

% Experimental glucose values for the GTT 
G_GTT_hyper=[10.096 9.997];
G_GTT_hyper=[10.22083333 11.17];
G_GTT_normo=[10.096 9.997];

index_G_hyper=[1 3];
index_G_normo=[2 4];
index_I_hyper=[5 7];
index_I_normo=[6 8];

% Experimental SEM values for the GTT
SEM_G_hyper=[0.652 0.652];
SEM_G_normo=[0.899 0.899];

SEM_G_hyper=[0 0];
SEM_G_normo=[0 0];

N_samples_G0=2;

% Index to select an insulin dose of 3500 mIU/L (included in the paper)
index=2;

% Select variables corresponding to the insulin dose
dose=dose(index);
G_GTT_hyper=G_GTT_hyper(index);
G_GTT_normo=G_GTT_normo(index);
index_G_hyper=index_G_hyper(index);
index_G_normo=index_G_normo(index);
index_I_hyper=index_I_hyper(index);
index_I_normo=index_I_normo(index);
SEM_G_hyper=SEM_G_hyper(index);
SEM_G_normo=SEM_G_normo(index);

% Define time axes to plot

simTime_5p5mM=[0:0.01:time_end];
simTime_11mM=[0:0.01:time_end];

nTime_11mM = length(simTime_11mM);

end_time = 336; % Plot for the total duration of co-culture time (336 h)

%% Hypothesis 1

Hypothesis1='H1_experiment';

optModel = IQMmodel(strcat(Hypothesis1,'.txt')); 
%IQMmakeMEXmodel(optModel,Hypothesis1); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);

loadParameterIndexesH1

file = ['parameterValues H1_experiment.dat']; 

tmpHold = load(file);

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    dataH1 = plotFunctionExperimentH1(X,Hypothesis1,parIndex,X_opt, simTime_11mM,...
        simTime_5p5mM,dose,limit_iiuptake,G_GTT_hyper,G_GTT_normo,N_samples_G0,SEM_G_hyper,SEM_G_normo);
end

%% Hypothesis 2

Hypothesis2='H2_experiment';

optModel = IQMmodel(strcat(Hypothesis2,'.txt')); 
%IQMmakeMEXmodel(optModel,Hypothesis2); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);

loadParameterIndexesH2

file = ['parameterValues H2_experiment.dat']; 

tmpHold = load(file);

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    dataH2 = plotFunctionExperimentH2(X,Hypothesis2,parIndex,X_opt, simTime_11mM,...
        simTime_5p5mM,dose,limit_iiuptake,ratio_F_hyper,G_GTT_hyper,G_GTT_normo);
end

title ('Glucose prediction with insulin dose for Hypothesis 1 (red) and Hypothesis 2 (green)')
%% Compare predictions of glucose and insulin with experimental data 
% for each hypothesis

% Plot predicted vs. experimental data. Predicted responses are represented
% as simulation lines (one line corresponds to an acceptable parameter set)

figure()

% Plot predicted vs. experimental data for hypothesis 1
% Plot predictions

% Plot hyperglycemic responses
for j = 1 : size(param,1)
    plot(simTime_11mM,dataH1.Gmeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    hold on
end

hold on

% Plot normoglycemic responses
for j = 1 : size(param,1)
    plot(simTime_11mM,dataH1.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM,'Linewidth',1.5)
    hold on
end

xlim([288 end_time])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose (mM)','FontSize',20)
xticks([0 48 96 144 192 240 288 336])
set(gca,'TickDir','out','FontSize',20);
box off

% Plot experimental data

% Plot hyperglycemic experimental data
hold on
errorbar(VALIDATIONDATA.time{index_G_hyper}, VALIDATIONDATA.mean{index_G_hyper}, VALIDATIONDATA.SD{index_G_hyper},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on

% Plot normoglycemic experimental data
hold on
errorbar(VALIDATIONDATA.time{index_G_normo}, VALIDATIONDATA.mean{index_G_normo}, VALIDATIONDATA.SD{index_G_normo},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on

title ('Hypothesis 1')

% Plot predicted vs. experimental data for hypothesis 2
figure()

% Plot predictions

% Plot hyperglycemic responses
for j = 1 : size(param,1)
    plot(simTime_11mM,dataH2.Gmeasured_11mM(j,1:nTime_11mM),'Color',color11mM,'Linewidth',1.5)
    hold on
end

hold on

% Plot normoglycemic responses
for j = 1 : size(param,1)
    plot(simTime_11mM,dataH2.Gmeasured_5p5mM(j,1:nTime_11mM),'Color',color5p5mM,'Linewidth',1.5)
    hold on
end

xlim([288 end_time])
ylim([0 13.75])
yticks([0 2.75 5.5 8.25 11 13.75])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
xticks([0 48 96 144 192 240 288 336])
box off

% Plot experimental data

% Plot hyperglycemic experimental data
hold on
errorbar(VALIDATIONDATA.time{index_G_hyper}, VALIDATIONDATA.mean{index_G_hyper}, VALIDATIONDATA.SD{index_G_hyper},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on

% Plot normoglycemic experimental data
hold on
errorbar(VALIDATIONDATA.time{index_G_normo}, VALIDATIONDATA.mean{index_G_normo}, VALIDATIONDATA.SD{index_G_normo},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on

title ('Hypothesis 2')

% Generate a plot with uncertainties represented as areas
% Calculate maximal and minimal values of the predictions to plot areas

% Hypothesis 1
dataH1.Gmeasured_11mM_max=max(dataH1.Gmeasured_11mM,[],1);
dataH1.Gmeasured_11mM_min=min(dataH1.Gmeasured_11mM,[],1);

dataH1.Gmeasured_5p5mM_max=max(dataH1.Gmeasured_5p5mM,[],1);
dataH1.Gmeasured_5p5mM_min=min(dataH1.Gmeasured_5p5mM,[],1);

dataH1.Imeasured_11mM_max=max(dataH1.Imeasured_11mM,[],1);
dataH1.Imeasured_11mM_min=min(dataH1.Imeasured_11mM,[],1);

dataH1.Imeasured_5p5mM_max=max(dataH1.Imeasured_5p5mM,[],1);
dataH1.Imeasured_5p5mM_min=min(dataH1.Imeasured_5p5mM,[],1);

dataH1.Gdifference=dataH1.Gmeasured_11mM-dataH1.Gmeasured_5p5mM;
dataH1.Idifference=dataH1.Gmeasured_11mM-dataH1.Gmeasured_5p5mM;

dataH1.Gdifference_max=max(dataH1.Gdifference,[],1);
dataH1.Gdifference_min=min(dataH1.Gdifference,[],1);

dataH1.Idifference_max=max(dataH1.Idifference,[],1);
dataH1.Idifference_min=min(dataH1.Idifference,[],1);

% Hypothesis 2
dataH2.Gmeasured_11mM_max=max(dataH2.Gmeasured_11mM,[],1);
dataH2.Gmeasured_11mM_min=min(dataH2.Gmeasured_11mM,[],1);

dataH2.Gmeasured_5p5mM_max=max(dataH2.Gmeasured_5p5mM,[],1);
dataH2.Gmeasured_5p5mM_min=min(dataH2.Gmeasured_5p5mM,[],1);

dataH2.Imeasured_11mM_max=max(dataH2.Imeasured_11mM,[],1);
dataH2.Imeasured_11mM_min=min(dataH2.Imeasured_11mM,[],1);

dataH2.Imeasured_5p5mM_max=max(dataH2.Imeasured_5p5mM,[],1);
dataH2.Imeasured_5p5mM_min=min(dataH2.Imeasured_5p5mM,[],1);

dataH2.Gdifference=dataH2.Gmeasured_11mM-dataH2.Gmeasured_5p5mM;
dataH2.Idifference=dataH2.Gmeasured_11mM-dataH2.Gmeasured_5p5mM;

dataH2.Gdifference_max=max(dataH2.Gdifference,[],1);
dataH2.Gdifference_min=min(dataH2.Gdifference,[],1);

dataH2.Idifference_max=max(dataH2.Idifference,[],1);
dataH2.Idifference_min=min(dataH2.Idifference,[],1);


%% Plot glucose predictions

% Plot difference in predictions between hypothesis 1 and 2
% figure()
% 
% for j = 1 : size(param,1)
%     plot(simTime_11mM,dataH1.Gdifference(j,1:nTime_11mM),'Color',color5p5mM,'Linewidth',1.5)
%     hold on
% end
% 
% title('Glucose difference')
% 
% % Calculate mean and SD of the difference between hyper and normoglycemia
% 
% mean_diff=VALIDATIONDATA.mean{index_G_hyper}-VALIDATIONDATA.mean{index_G_normo};
% SD_diff=sqrt(VALIDATIONDATA.SD{index_G_hyper}.^2/...
%     N_replicates+VALIDATIONDATA.SD{index_G_normo}.^2/N_replicates);
% 
% hold on
% errorbar(VALIDATIONDATA.time{index_G_hyper}, ...
%     mean_diff,...
%     SD_diff,'Color', color11mM,...
%     'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
% 
% xlim([288 end_time])

% Plot glucose and insulin responses in comparison to experimental data

figure()

hold on

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[dataH1.Gmeasured_5p5mM_max,fliplr(dataH1.Gmeasured_5p5mM_min)];
f=fill(area1,area2,color5p5mM,'LineStyle','none');
alpha(f,.2)

hold on

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[dataH1.Gmeasured_11mM_max,fliplr(dataH1.Gmeasured_11mM_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(VALIDATIONDATA.time{index_G_hyper}, ...
    VALIDATIONDATA.mean{index_G_hyper}, VALIDATIONDATA.SD{index_G_hyper},...
    'Color', color11mM,'LineStyle', 'n', 'Marker',marker,...
    'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on

% Glucose 5.5 mM, dose A
hold on
errorbar(VALIDATIONDATA.time{index_G_normo}, ...
    VALIDATIONDATA.mean{index_G_normo}, VALIDATIONDATA.SD{index_G_normo},...
    'Color', color5p5mM,'LineStyle', 'n', 'Marker',marker,...
    'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on


xlim([288 end_time])
ylim([0 13.75])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
xticks([0 48 96 144 192 240 288 336])
yticks([0 2.75 5.5 8.25 11 13.75])
box off
title('Hypothesis 1')


% Plot prediction for hypothesis 2

figure()

hold on

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[dataH2.Gmeasured_5p5mM_max,fliplr(dataH2.Gmeasured_5p5mM_min)];
f=fill(area1,area2,color5p5mM,'LineStyle','none');
alpha(f,.2)

hold on

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[dataH2.Gmeasured_11mM_max,fliplr(dataH2.Gmeasured_11mM_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(VALIDATIONDATA.time{index_G_hyper}, ...
    VALIDATIONDATA.mean{index_G_hyper}, ...
    VALIDATIONDATA.SD{index_G_hyper},'Color', color11mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on

% Glucose 5.5 mM, dose A
hold on
errorbar(VALIDATIONDATA.time{index_G_normo}, ...
    VALIDATIONDATA.mean{index_G_normo}, ...
    VALIDATIONDATA.SD{index_G_normo},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on


xlim([288 end_time])
ylim([0 13.75])
xlabel('Time (h)','FontSize',20);
ylabel('Glucose (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
xticks([0 48 96 144 192 240 288 336])
yticks([0 2.75 5.5 8.25 11 13.75])
box off
title('Hypothesis 2')

%% Plot insulin predictions

% Plot prediction for hypothesis 1

figure()

hold on

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[dataH1.Imeasured_5p5mM_max,fliplr(dataH1.Imeasured_5p5mM_min)];
f=fill(area1,area2,color5p5mM,'LineStyle','none');
alpha(f,.2)

hold on

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[dataH1.Imeasured_11mM_max,fliplr(dataH1.Imeasured_11mM_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(VALIDATIONDATA.time{index_I_hyper}, ...
    VALIDATIONDATA.mean{index_I_hyper}, VALIDATIONDATA.SD{index_I_hyper},...
    'Color', color11mM,'LineStyle', 'n', 'Marker',marker,...
    'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on

% Glucose 5.5 mM, dose A
hold on
errorbar(VALIDATIONDATA.time{index_I_normo}, ...
    VALIDATIONDATA.mean{index_I_normo}, ...
    VALIDATIONDATA.SD{index_I_normo},'Color', color5p5mM,...
    'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on


xlim([288 simTime_11mM(end)])
ylim([0 4000])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
xticks([0 48 96 144 192 240 288 336])
box off
title('Hypothesis 1')

% Plot prediction for hypothesis 2

figure()

hold on

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[dataH2.Imeasured_5p5mM_max,fliplr(dataH2.Imeasured_5p5mM_min)];
f=fill(area1,area2,color5p5mM,'LineStyle','none');
alpha(f,.2)

hold on

area1=[simTime_11mM,fliplr(simTime_11mM)];
area2=[dataH2.Imeasured_11mM_max,fliplr(dataH2.Imeasured_11mM_min)];
f=fill(area1,area2,color11mM,'LineStyle','none');
alpha(f,.2)

hold on
errorbar(VALIDATIONDATA.time{index_I_hyper}, ...
    VALIDATIONDATA.mean{index_I_hyper}, VALIDATIONDATA.SD{index_I_hyper},...
    'Color', color11mM,'LineStyle', 'n', 'Marker',marker,...
    'MarkerFaceColor',color11mM, 'MarkerSize',8)
hold on

% Glucose 5.5 mM, dose A
hold on
errorbar(VALIDATIONDATA.time{index_I_normo}, ...
    VALIDATIONDATA.mean{index_I_normo}, VALIDATIONDATA.SD{index_I_normo},...
    'Color', color5p5mM,'LineStyle', 'n', 'Marker',marker,...
    'MarkerFaceColor',color5p5mM, 'MarkerSize',8)
hold on


xlim([288 simTime_11mM(end)])
ylim([0 4000])
xlabel('Time (h)','FontSize',20);
ylabel('Insulin (mM)','FontSize',20)
set(gca,'TickDir','out','FontSize',20);
xticks([0 48 96 144 192 240 288 336])
box off
title('Hypothesis 2')

[time_start index_start]=find(simTime_11mM==288);

simTime_11mM=simTime_11mM(index_start:end);

% Hypothesis 1

dataH1.Gmeasured_11mM_max=dataH1.Gmeasured_11mM_max(index_start:end)';
dataH1.Gmeasured_11mM_min=dataH1.Gmeasured_11mM_min(index_start:end)';

dataH1.Gmeasured_5p5mM_max=dataH1.Gmeasured_5p5mM_max(index_start:end)';
dataH1.Gmeasured_5p5mM_min=dataH1.Gmeasured_5p5mM_min(index_start:end)';

dataH1.Imeasured_11mM_max=dataH1.Imeasured_11mM_max(index_start:end)';
dataH1.Imeasured_11mM_min=dataH1.Imeasured_11mM_min(index_start:end)';

dataH1.Imeasured_5p5mM_max=dataH1.Imeasured_5p5mM_max(index_start:end)';
dataH1.Imeasured_5p5mM_min=dataH1.Imeasured_5p5mM_min(index_start:end)';

dataH1.Gdifference=(dataH1.Gmeasured_11mM(:,index_start:end)-...
    dataH1.Gmeasured_5p5mM(:,index_start:end));
dataH1.Idifference=(dataH1.Gmeasured_11mM(:,index_start:end)-...
    dataH1.Gmeasured_5p5mM(:,index_start:end));

dataH1.Gdifference_max=max(dataH1.Gdifference,[],1)';
dataH1.Gdifference_min=min(dataH1.Gdifference,[],1)';

dataH1.Idifference_max=max(dataH1.Idifference,[],1)';
dataH1.Idifference_min=min(dataH1.Idifference,[],1)';

% Hypothesis 2

dataH2.Gmeasured_11mM_max=dataH2.Gmeasured_11mM_max(index_start:end)';
dataH2.Gmeasured_11mM_min=dataH2.Gmeasured_11mM_min(index_start:end)';

dataH2.Gmeasured_5p5mM_max=dataH2.Gmeasured_5p5mM_max(index_start:end)';
dataH2.Gmeasured_5p5mM_min=dataH2.Gmeasured_5p5mM_min(index_start:end)';

dataH2.Imeasured_11mM_max=dataH2.Imeasured_11mM_max(index_start:end)';
dataH2.Imeasured_11mM_min=dataH2.Imeasured_11mM_min(index_start:end)';

dataH2.Imeasured_5p5mM_max=dataH2.Imeasured_5p5mM_max(index_start:end)';
dataH2.Imeasured_5p5mM_min=dataH2.Imeasured_5p5mM_min(index_start:end)';

dataH2.Gdifference=(dataH2.Gmeasured_11mM(:,index_start:end)-...
    dataH2.Gmeasured_5p5mM(:,index_start:end));
dataH2.Idifference=(dataH2.Gmeasured_11mM(:,index_start:end)-...
    dataH2.Gmeasured_5p5mM(:,index_start:end));

dataH2.Gdifference_max=max(dataH2.Gdifference,[],1)';
dataH2.Gdifference_min=min(dataH2.Gdifference,[],1)';

dataH2.Idifference_max=max(dataH2.Idifference,[],1)';
dataH2.Idifference_min=min(dataH2.Idifference,[],1)';

% Compute statistical chi-2 tests for rejecting/accepting hypothesis 

% Number of data points in the validation dataset. The first sample in the
% GTT (glucose and insulin at t=0) is not included in the calculation

N=length(VALIDATIONDATA.time{index_G_hyper}(2:end))...
    +length(VALIDATIONDATA.time{index_G_normo}(2:end));

chi2_threshold=chi2inv(0.95,N);

errorH1=zeros(1,size(param,1));
simTime_11mM=[0:0.01:time_end];
% Find time values in the validation dataset within the simulated times
time_index=find(ismember(simTime_11mM,VALIDATIONDATA.time{index_G_hyper}));

% Compute error for each parameter set for H1
for j = 1 : size(param,1)

    errorH1(j)=sum(((dataH1.Gmeasured_11mM(j,time_index)'- ...
        VALIDATIONDATA.mean{index_G_hyper}).^2)./(VALIDATIONDATA.SD{index_G_hyper}).^2)+...
        sum(((dataH1.Gmeasured_5p5mM(j,time_index)'- ...
        VALIDATIONDATA.mean{index_G_normo}).^2)./(VALIDATIONDATA.SD{index_G_normo}).^2);
end

% Minimal error accross all parameter values, to compare with the chi-2
% value
minerrorH1=min(errorH1);

% Compute the error for each parameter set for H2
errorH2=zeros(1,size(param,1));

% Compute the error for each parameter set for H1
for j = 1 : size(param,1)

    errorH2(j)=sum(((dataH2.Gmeasured_11mM(j,time_index)'- ...
        VALIDATIONDATA.mean{index_G_hyper}).^2)./...
        (VALIDATIONDATA.SD{index_G_hyper}).^2)+...
        sum(((dataH2.Gmeasured_5p5mM(j,time_index)'- ...
        VALIDATIONDATA.mean{index_G_normo}).^2)./(VALIDATIONDATA.SD{index_G_normo}).^2);
end

minerrorH2=min(errorH2);
