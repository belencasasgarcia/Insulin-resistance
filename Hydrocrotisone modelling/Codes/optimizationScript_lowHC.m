% Optimization script: 
% 
% This script estimates the set of parameters that
% provide an acceptable agreement between the simulated glucose and insulin
% responses from the model corresponding to Hypothesis 2, considering low 
% hydrocortisone (HC), during the first 48 h of co-culture, 
% and the corresponding experimental data. These parameter values are used
% for baseline correction of the predictions at low HC.


clear
clc
close all

format long
format compact

% Declare global variables

global EXPDATA
global modelName
global icOrig0
global pNames
global param
global OPTIONS
global COSTOPTIONS
global SIMTIME      % Time axis for simulation
global time_highres % Time axis for plotting simulations 
global OPTVAR_INDEX    % Indexes of variables to be optimized
global CUTOFF
global dataPoints
global FID
global Optparam
global parIndex


%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData_HC

modelName ='H2_experiment';

optModel = IQMmodel(strcat(modelName,'.txt')); 
%IQMmakeMEXmodel(optModel,modelName); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);

% figure()
% plotSimulation(param)

%% Parameter indexes

loadParameterIndexes


%% Optimization settings

% Restrict optimization to the first 48 hours
time_end=48;

SIMTIME=[0:0.1:time_end];

time_highres=[0:0.001:time_end];

% Find experimental values corresponding to the simulation time

for i=1:size(EXPDATA.time,2)
    time_indexes=find(EXPDATA.time{i}<=time_end);
    EXPDATA.time{i}=EXPDATA.time{i}(time_indexes);
    EXPDATA.mean{i}=EXPDATA.mean{i}(time_indexes);
    EXPDATA.SD{i}=EXPDATA.SD{i}(time_indexes);
end

N_iter=10;             % Number of iterations for the optimization

Optparam=param;        %Optimal parameters

OPTVAR_INDEX=[7 8];    %Experimental data to optimize. 
OPTIONS.lowbounds=[];
OPTIONS.highbounds=[];

% Define parameters to optimize. The parameters adjusted for baseline
% correction are insulin sensitivity (S_i), insulin secretion capacity (Sigma)
% and glucose and insulin offsets at t=0 (start of the GTT)

OPTIONS.lowbounds(1)=1e-9;                % S_i
OPTIONS.lowbounds(2)=1e6;                 % Sigma
OPTIONS.lowbounds(3)=-2;                  % delta_G_1_11mM
OPTIONS.lowbounds(4)=0;                   % delta_I_1_11mM


OPTIONS.lowbounds=OPTIONS.lowbounds';

% High parameter bounds

OPTIONS.highbounds(1)=1;                  % S_i
OPTIONS.highbounds(2)=1e9;                % Sigma
OPTIONS.highbounds(3)=2;                  % delta_G_1_11mM
OPTIONS.highbounds(4)=150;                % delta_I_1_11mM

OPTIONS.highbounds=OPTIONS.highbounds';

% Get startGuess for the parameters from the optimization at high HC


startGuess(1)=0.003;                      % S_i
startGuess(2)=13e6;                       % Sigma
startGuess(3)=0;                          % delta_G_1_11mM
startGuess(4)=0;                          % delta_I_1_11mM


OPTIONS.index_Optpar=logical(sum([parIndex.i_S_i parIndex.i_Sigma ...
    parIndex.i_delta_G_1_11mM ...
    parIndex.i_delta_I_1_11mM],2)); %Indexes of the parameters to optimize

% The start guess for the parameters are the optimal parameter values from
% high HCT. These parameters are saved as 'popt_highHC.mat'

popt = load('popt_H2_experiment_highHC.mat');
Optparam = popt.X_opt;

startGuess = Optparam(OPTIONS.index_Optpar)

X=startGuess;

%Set parameter values manually

Optparam(parIndex.i_EGP_hep)=0;
Optparam(parIndex.i_G0)=11;
Optparam(parIndex.i_G_GTT)=11;

icOrig0=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
   0*param(parIndex.i_V_m_hep) 0*param(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];

% Select start guess as initial value for the parameters

Optparam(OPTIONS.index_Optpar)=startGuess;

startCost = costFunction_lowHC(startGuess);

% Simulated annealing
% Set temperature
OPTIONS.tempstart = 1e1*startCost;                   % InitialTemp
OPTIONS.maxitertemp =   50*length(startGuess);       % Max-iterations per temp
OPTIONS.maxitertemp0 =  200*length(startGuess);      % Max-iterations per temp0c

% SBPDsimulate optionsc
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through

% Costfunction settings
COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

%% Choose data to optimize

%  Cut-off value for chi-2 test

% Calculate number of data points for statistical tests
dataPoints=0;

for i=1:OPTVAR_INDEX(end)
    dataPoints=dataPoints+length(EXPDATA.time{i});
end

CUTOFF = chi2inv(0.99,dataPoints);

file = ['parameterValues' ' ' modelName ' ' datestr(datetime('now')) '.dat'];

FID = fopen(file, 'wt');

N_iter=1;      % Number of iterations for the optimization

for i=1:N_iter
    startGuess=X;
    startCost = costFunction_lowHC(startGuess);
    disp(['Initial cost with start guess: ' num2str(startCost)]);
    
    [X,FVAL,EXITFLAG] = simannealingSBAOClusteringL(@costFunction_lowHC,startGuess,OPTIONS);
    
end


