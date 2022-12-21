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

modelName ='M1_experiment';

optModel = IQMmodel(strcat(modelName,'.txt')); 
IQMmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

time_end=336;

SIMTIME=[0:0.005:time_end];

% Parameter indexes

%% Parameter indexes

loadParameterIndexes

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


file = ['parameterValues M1_experiment 30-Jul-2021 12:03:27.dat'];
%% 

tmpHold = load(file);

% Define time axes to plot

simTime_5p5mM=[0:0.01:time_end];
simTime_11mM=[0:0.01:time_end];

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    plotFunctionExperiment(X,modelName,parIndex,X_opt, simTime_11mM,...
        simTime_5p5mM)
end




