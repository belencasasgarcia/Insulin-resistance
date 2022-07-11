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

modelName ='Islets_liver_2OC';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

SIMTIME=[0:0.001:48];

% Parameter indexes

parIndex.i_V_m_hep = ismember(pNames,'V_m_hep');
parIndex.i_V_hep = ismember(pNames,'V_hep');
parIndex.i_V_m_islets = ismember(pNames,'V_m_islets');
parIndex.i_V_islets = ismember(pNames,'V_islets');
parIndex.i_Q = ismember(pNames,'Q');
parIndex.i_EGP_hep = ismember(pNames,'EGP_hep');
parIndex.i_U_ii_hep = ismember(pNames,'U_ii_hep');
parIndex.i_S_i = ismember(pNames,'S_i');
parIndex.i_Sigma = ismember(pNames,'Sigma');
parIndex.i_Alpha = ismember(pNames,'Alpha');
parIndex.i_CL_i_hep = ismember(pNames,'CL_i_hep');
parIndex.i_V_sample_hep = ismember(pNames,'V_sample_hep');
parIndex.i_V_sample_islets = ismember(pNames,'V_sample_islets')

iC=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
    183*param(parIndex.i_V_m_hep) 183*param(parIndex.i_V_m_islets)];

%iC=[2.8*param(parIndex.i_V_m_hep) 2.8*param(parIndex.i_V_m_islets)...
%    0*param(parIndex.i_V_m_hep) 0*param(parIndex.i_V_m_islets)];

% SBPDsimulate options
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through

file = ['parameterValues Islets_liver_2OC 21-Oct-2019 10:57:14.dat'];

tmpHold = load(file);


if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    plotFunction(X,'Islets_liver_2OC',SIMTIME,iC)
end
