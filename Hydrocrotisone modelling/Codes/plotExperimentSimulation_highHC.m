% Plot model simulations with uncertainty based on the estimated parameter
% values

clear
clc
close all
format long
format compact

% Declare global variables

global modelName
global icOrig0
global pNames
global param
global OPTIONS
global SIMTIME      % Time axis for simulation
global parIndex

%% Load parameter indexes

loadParameterIndexes

%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData_HC

modelName ='H2_experiment';

optModel = IQMmodel(strcat(modelName,'.txt')); 
%SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);

SIMTIME=[0:0.005:336];

% Set initial conditions for the simulation

iC=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
   0*param(parIndex.i_V_m_hep) 0*param(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088];

% SBPDsimulate options
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through

% Load acceptable parameters from the optimization routine

file = ['parameterValues H2_experiment_highHC.dat'];

tmpHold = load(file);

% Define time axes to plot

% Time axis for normoglycemia
simTime_5p5mM=[0:0.05:336];

% Time axis for hyperglycemia
simTime_11mM=[0:0.05:336];

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    plotFunctionUncertainty(X,modelName,parIndex,X_opt, simTime_11mM,...
        simTime_5p5mM)
end




