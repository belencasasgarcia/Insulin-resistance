% Plot model simulations with uncertainty based on the estimated parameter
% values

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
global OPTIME       % Time axis for optimization
global N_VAR        % Number of variables to optimize
global VAR_INDEX    % Indexes of variables to be optimized
global CUTOFF
global dataPoints
global FID
global Optparam
global parIndex

%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData

modelName ='M3';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

SIMTIME=[0:0.005:336];

% Parameter indexes

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
parIndex.i_k = ismember(pNames,'k');
parIndex.i_kv = ismember(pNames,'kv');
parIndex.i_Vm_f_add = ismember(pNames,'Vm_f_add');
parIndex.i_Km_f_add = ismember(pNames,'Km_f_add');
parIndex.i_Vm_f_islets = ismember(pNames,'Vm_f_islets');
parIndex.i_Km_f_islets = ismember(pNames,'Km_f_islets');

% Media change parameters
parIndex.i_G0 = ismember(pNames,'G0');
parIndex.i_I0 = ismember(pNames,'I0');
parIndex.i_G_GTT = ismember(pNames,'G_GTT');
parIndex.i_I_GTT = ismember(pNames,'I_GTT');
parIndex.i_G0_7 = ismember(pNames,'G0_7');
parIndex.i_I0_7 = ismember(pNames,'I0_7');

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


file = ['parameterValues M1 26-Jul-2021 18:22:45.dat'];

tmpHold = load(file);

% Define time axes to plot

simTime_5p5mM=[0:0.01:192];
simTime_11mM=[0:0.01:192];

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    plotFunctionUncertainty_lines(X,modelName,parIndex,X_opt, simTime_11mM,...
        simTime_5p5mM)
end




