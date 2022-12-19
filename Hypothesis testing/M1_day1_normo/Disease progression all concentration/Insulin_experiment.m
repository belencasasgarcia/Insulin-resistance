% Add insulin

clear
clc
close all

format long
format compact

% Declare global variables

global EXPDATA_G
global EXPDATA_I
global EXPDATA_G_hep
global modelName
global icOrig0
global pNames
global param
global OPTIONS
global COSTOPTIONS
global EXPDATA
global SIMTIME      % Time axis for simulation
global time_highres % Time axis for plotting simulations 
global OPTIME       % Time axis for optimization
global N_VAR        % Number of variables to optimize
global OPTVAR_INDEX    % Indexes of variables to be optimized
global CUTOFF
global dataPoints
global FID
global Optparam
global parIndex


%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData

modelName ='Islets_liver_2OC_dis_prog';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

% figure()
% plotSimulation(param)

%% Parameter indexes

parIndex.i_V_m_hep = ismember(pNames,'V_m_hep');
parIndex.i_V_hep = ismember(pNames,'V_hep');
parIndex.i_V_m_islets = ismember(pNames,'V_m_islets');
parIndex.i_Q = ismember(pNames,'Q');
parIndex.i_EGP_hep = ismember(pNames,'EGP_hep');
parIndex.i_U_ii_hep = ismember(pNames,'U_ii_hep');
parIndex.i_S_i = ismember(pNames,'S_i');
parIndex.i_Sigma = ismember(pNames,'Sigma');
parIndex.i_Alpha = ismember(pNames,'Alpha');
parIndex.i_CL_i_hep = ismember(pNames,'CL_i_hep');
parIndex.i_V_sample_hep = ismember(pNames,'V_sample_hep');
parIndex.i_V_sample_islets = ismember(pNames,'V_sample_islets');

% Disease progession parameters
parIndex.i_Vm = ismember(pNames,'Vm');
parIndex.i_km = ismember(pNames,'km');
parIndex.i_d0 = ismember(pNames,'d0');
parIndex.i_r1 = ismember(pNames,'r1');
parIndex.i_r2 = ismember(pNames,'r2');
parIndex.i_tao_slow = ismember(pNames,'tao_slow');
parIndex.i_G_healthy = ismember(pNames,'G_healthy');
parIndex.i_hep = ismember(pNames,'hep');
parIndex.i_islets = ismember(pNames,'islets');

% Media change parameters
parIndex.i_G0 = ismember(pNames,'G0');
parIndex.i_I0 = ismember(pNames,'I0');
parIndex.i_G_GTT = ismember(pNames,'G_GTT');

parIndex.i_k = ismember(pNames,'k');
%% Optimization settings

% Simulation time
time_end=336;

SIMTIME=[0:0.01:time_end];
%SIMTIME=[0:0.01:48];

time_highres=[0:0.001:time_end];

N_iter=1;      % Number of iterations for the optimization

Optparam=param;  %Optimal parameters

parameters=load('Optpar_allconc.mat')

Optparam=parameters.Optparam;
 
OPTIONS.index_Optpar=boolean(sum([parIndex.i_U_ii_hep...
    parIndex.i_Sigma parIndex.i_CL_i_hep parIndex.i_Vm...
    parIndex.i_km parIndex.i_d0 parIndex.i_r1 parIndex.i_r2 ...
    parIndex.i_k],2)); %Indexes of the parameters to optimize

%Set parameter values manually
Optparam(parIndex.i_EGP_hep)=0;
Optparam(parIndex.i_G0)=11;
Optparam(parIndex.i_G_GTT)=11;
Optparam(parIndex.i_U_ii_hep)=0.45;

icOrig0=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
   183.3*param(parIndex.i_V_m_hep) 183.3*param(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088];

% Select start guess as initial value for the parameters

% Plot simulation and data corresponding to the initial guess
SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig0,pNames,Optparam);

figure()
plotSimulation(SSsimData,0)

% Plot function defining beta cell growth
plot_beta_dynamics(Optparam(parIndex.i_d0),Optparam(parIndex.i_r1),...
    Optparam(parIndex.i_r2))


%% Plot the result of adding insulin (no islets are present)
% One has to include different initial conditions depending on the specific
% concentration

icOrig0_exp=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
   800*param(parIndex.i_V_m_islets) 800*param(parIndex.i_V_m_islets) 0 9.65e+02 ...
   6.98 0.0000000088];

modelName_exp='Islets_liver_2OC_dis_prog_insulin';
SIMTIME_exp=[0:0.01:48];

expModel = SBmodel(strcat(modelName_exp,'.txt')); 
SBPDmakeMEXmodel(expModel,modelName_exp); 

Optparam(parIndex.i_islets)=0;
Optparam(parIndex.i_CL_i_hep)=2;
Optparam(parIndex.i_U_ii_hep)=1;
Optparam(parIndex.i_S_i)=0.006593714275828;

% Plot simulation
SimData_exp=SBPDsimulate(modelName_exp,SIMTIME_exp,icOrig0_exp,pNames,Optparam);

figure()
plotSimulation_exp(SimData_exp,0,SIMTIME_exp,expModel)








