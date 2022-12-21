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

modelName ='M3_experiment';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

SIMTIME=[0:0.005:336];

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


%file = ['parameterValues M3_experiment 02-Aug-2021 15:50:25.dat'];
%file = ['parameterValues M3_experiment 02-Aug-2021 13:28:53.dat'];
%file = ['parameterValues M3_experiment 31-Jul-2021 17:32:25.dat'];parameterValues M3_experiment 02-Aug-2021 15:43:04.dat
%file = ['parameterValues M3_experiment 31-Jul-2021 00:35:40.dat'];

file = ['parameterValues M3_experiment 30-Aug-2021 17:10:15.dat']; %Works!!!

%file = ['parameterValues M3_experiment 08-Sep-2021 10:14:51.dat']; %Works!!!
%% 

tmpHold = load(file);

% Define time axes to plot

simTime_5p5mM=[0:0.01:191.9];
simTime_11mM=[0:0.01:191.9];

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    plotFunctionExperiment(X,modelName,parIndex,X_opt, simTime_11mM,...
        simTime_5p5mM)
end




