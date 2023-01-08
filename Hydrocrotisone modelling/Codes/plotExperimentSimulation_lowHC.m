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

%% Load parameter indexes

loadParameterIndexes

%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load experimental data

loadExperimentalData_HC

modelName ='H2_experiment_lowHC';

optModel = IQMmodel(strcat(modelName,'.txt')); 
%IQMmakeMEXmodel(optModel,modelName); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);

SIMTIME=[0:0.005:336];

% Parameter indexes
loadParameterIndexes

iC=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
   183.3*param(parIndex.i_V_m_hep) 183.3*param(parIndex.i_V_m_islets) 0 0 ...
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

file = ['parameterValues H2_experiment_highHC.dat'];

tmpHold = load(file);

% Define time axes to plot

simTime_5p5mM=[0:0.05:360];

simTime_11mM=[0:0.05:360];

% Load datasets for the parameters between days 1 and 3

indexes_optimized_GTTd1to3_spreaded = load('indexes_optimized_GTTd1to3_spreaded.mat');
indexes_optimized_GTTd1to3_spreaded = indexes_optimized_GTTd1to3_spreaded.indexes_optimized_GTTd1to3_spreaded;

opt_par_GTTd1to3_spreaded = load('opt_par_GTTd1to3_spreaded.mat');
opt_par_GTTd1to3_spreaded = opt_par_GTTd1to3_spreaded.X;

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    % Change offsets
    X_opt = tmpHold(end,2:end);
    plotFunctionUncertainty_lowHC(X,modelName,parIndex,X_opt, simTime_11mM,...
        simTime_5p5mM)
end





