% Main simulation file 

% You have changed the range for km in the insulin resistance parameters
% and introduced sample 2 in the spheroid data.

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

loadExperimentalDataExp1_corrected_lowHC

modelName ='M3_experiment';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

% figure()
% plotSimulation(param)

%% Parameter indexes

loadParameterIndexes


%% Optimization settings

% Simulation time
time_end=48;

SIMTIME=[0:0.1:time_end];
%SIMTIME=[0:0.01:48];

time_highres=[0:0.001:time_end];

% Change experimental data based on the simulation time

% Find experimental values corresponding to the simulation time

for i=1:size(EXPDATA.time,2)
    time_indexes=find(EXPDATA.time{i}<=time_end);
    EXPDATA.time{i}=EXPDATA.time{i}(time_indexes);
    EXPDATA.mean{i}=EXPDATA.mean{i}(time_indexes);
    EXPDATA.SD{i}=EXPDATA.SD{i}(time_indexes);
end

plotExperimentalData_lowHC(EXPDATA)

N_iter=1;      % Number of iterations for the optimization

Optparam=param;  %Optimal parameters

OPTVAR_INDEX=[7 8];    %Experimental data to optimize. Keep this before settings!
OPTIONS.lowbounds=[];
OPTIONS.highbounds=[];

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


 
OPTIONS.index_Optpar=boolean(sum([parIndex.i_S_i parIndex.i_Sigma ...
    parIndex.i_delta_G_1_11mM ...
    parIndex.i_delta_I_1_11mM],2)); %Indexes of the parameters to optimize

popt = load('popt_highHC.mat');
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

% Plot simulation and data corresponding to the initial guess
SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig0,pNames,Optparam);

% Plot simulation with initial parameter values
figure()
plotSimulation_lowHC(SSsimData,0)

startCost = costFunction_lowHC(startGuess);

psOpts = optimoptions(@particleswarm,'Display','iter');

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

% Optimal parameter vector

Optparam(OPTIONS.index_Optpar)=X(1,:);

SimDataOpt=SBPDsimulate(modelName,SIMTIME,icOrig0,pNames,Optparam);

figure()
plotSimulation(SimDataOpt,0)

