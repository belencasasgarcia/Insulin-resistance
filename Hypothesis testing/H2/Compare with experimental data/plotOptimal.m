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

modelName ='Islets_liver_2OC_dis_prog';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

SIMTIME=[0:0.1:336];

iC=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
    183*param(parIndex.i_V_m_hep) 183*param(parIndex.i_V_m_islets) 0 1 ...
    0 5.5 0.0000000088];

% SBPDsimulate options
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through

optimal= load('Goodpar.mat');
parameters= optimal.Optparam;

plotSimulation(parameters,iC,0)

%% Simulate for 5.5 mM

iC=[5.5*param(parIndex.i_V_m_hep) 5.5*param(parIndex.i_V_m_islets)...
    0*param(parIndex.i_V_m_hep) 0*param(parIndex.i_V_m_islets) 0 1 ...
    0 5.5 0.0000000088];

parameters(parIndex.i_G0) = 5.5;

plotSimulation(parameters,iC,0)





