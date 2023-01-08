% Optimization script

% This script optimizes the model parameters against the experimental data.
% The result is a set of acceptable parameters, that provide an acceptable
% agreement with the experimental data according to a statistical chi-2
% test.

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

loadExperimentalData

modelName ='H2_experiment';

optModel = IQMmodel(strcat(modelName,'.txt')); 
% Uncomment the following line to create MEX model
%IQMmakeMEXmodel(optModel,modelName); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);

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
parIndex.i_I_GTT_13 = ismember(pNames,'I_GTT_13');

parIndex.i_delta_G_1 = ismember(pNames,'delta_G_1');
parIndex.i_delta_G_7 = ismember(pNames,'delta_G_7');

parIndex.i_delta_I_1 = ismember(pNames,'delta_I_1');
parIndex.i_delta_I_7 = ismember(pNames,'delta_I_7');

% Media change parameters for each concentration

% 11 mM

parIndex.i_delta_G_1_11mM = ismember(pNames,'delta_G_1_11mM');
parIndex.i_delta_G_7_11mM = ismember(pNames,'delta_G_7_11mM');

parIndex.i_delta_I_1_11mM = ismember(pNames,'delta_I_1_11mM');
parIndex.i_delta_I_7_11mM = ismember(pNames,'delta_I_7_11mM');

% 5.5 mM

parIndex.i_delta_G_7_5p5mM = ismember(pNames,'delta_G_7_5p5mM');
parIndex.i_delta_I_7_5p5mM = ismember(pNames,'delta_I_7_5p5mM');


%% Optimization settings

% Simulation time
time_end=288;

SIMTIME=[0:0.1:time_end];

time_highres=[0:0.001:time_end];

% Change experimental data based on the simulation time

% Find experimental values corresponding to the simulation time

for i=1:size(EXPDATA.time,2)
    time_indexes=find(EXPDATA.time{i}<=time_end);
    EXPDATA.time{i}=EXPDATA.time{i}(time_indexes);
    EXPDATA.mean{i}=EXPDATA.mean{i}(time_indexes);
    EXPDATA.SD{i}=EXPDATA.SD{i}(time_indexes);
end

N_iter=10;      % Number of iterations for the optimization

Optparam=param;  %Optimal parameters

OPTVAR_INDEX=[1 2 5 6];    %Experimental data to optimize. Glucose and 
% insulin for normo and hyperglycemia

OPTIONS.lowbounds=[];
OPTIONS.highbounds=[];

OPTIONS.lowbounds(1)=1e-2;                % U_ii
OPTIONS.lowbounds(2)=1e-9;                % S_i
OPTIONS.lowbounds(3)=1e6;                 % Sigma
OPTIONS.lowbounds(4)=1e-9;                % CL_i_hep
OPTIONS.lowbounds(5)=1e-9;                % Vm
OPTIONS.lowbounds(6)=200;                 % km
OPTIONS.lowbounds(7)=1e-9;                % kv
OPTIONS.lowbounds(8)=100;                 % Km_f_add
OPTIONS.lowbounds(9)=1;                   % Km_f_islets
OPTIONS.lowbounds(10)=-1.5;               % delta_G_1_11mM
OPTIONS.lowbounds(11)=-1.5;               % delta_G_7_11mM
OPTIONS.lowbounds(12)=0;                  % delta_I_1_11mM
OPTIONS.lowbounds(13)=0;                  % delta_I_7_11mM


OPTIONS.lowbounds=OPTIONS.lowbounds';

% High parameter bounds

OPTIONS.highbounds(1)=0.9;                % U_ii_hep
OPTIONS.highbounds(2)=0.1;                % S_i
OPTIONS.highbounds(3)=1e9;                % Sigma
OPTIONS.highbounds(4)=200;                % CL_i_hep
OPTIONS.highbounds(5)=1;                  % Vm
OPTIONS.highbounds(6)=1e5;                % km
OPTIONS.highbounds(7)=100;                % kv
OPTIONS.highbounds(8)=1500;               % Km_f_add
OPTIONS.highbounds(9)=500;                % Km_f_islets
OPTIONS.highbounds(10)=1.5;               % delta_G_1_11mM
OPTIONS.highbounds(11)=1.5;               % delta_G_7_11mM
OPTIONS.highbounds(12)=150;               % delta_I_1_11mM
OPTIONS.highbounds(13)=150;               % delta_I_7_11mM

OPTIONS.highbounds=OPTIONS.highbounds';

startGuess(1)=0.1;                        % U_ii_hep
startGuess(2)=0.003;                      % S_i
startGuess(3)=13e6;                       % Sigma
startGuess(4)=20;                         % CL_i_hep 
startGuess(5)=0.92;                       % Vm  
startGuess(6)=1200;                       % km 
startGuess(7)=0.001;                      % kv
startGuess(8)=200;                        % Km_f_add
startGuess(9)=50;                         % Km_f_islets
startGuess(10)=0;                         % delta_G_1_11mM
startGuess(11)=0;                         % delta_G_7_11mM
startGuess(12)=0;                         % delta_I_1_11mM
startGuess(13)=0;                         % delta_I_7_11mM


X=startGuess;
 
OPTIONS.index_Optpar=logical(sum([parIndex.i_U_ii_hep...
    parIndex.i_S_i parIndex.i_Sigma parIndex.i_CL_i_hep...
    parIndex.i_Vm parIndex.i_km...
    parIndex.i_kv parIndex.i_Km_f_add parIndex.i_Km_f_islets ...
    parIndex.i_delta_G_1_11mM parIndex.i_delta_G_7_11mM...
    parIndex.i_delta_I_1_11mM parIndex.i_delta_I_7_11mM],2)); %Indexes of the parameters to optimize

%Set parameter values manually
Optparam(parIndex.i_EGP_hep)=0;
Optparam(parIndex.i_G0)=11;
Optparam(parIndex.i_G_GTT)=11;

icOrig0=[9.519*param(parIndex.i_V_m_hep) 9.519*param(parIndex.i_V_m_islets)...
   131.1*param(parIndex.i_V_m_hep) 131.1*param(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];

% Select start guess as initial value for the parameters

Optparam(OPTIONS.index_Optpar)=startGuess;

startCost = costFunction(startGuess);

psOpts = optimoptions(@particleswarm,'Display','iter');

% Simulated annealing
% Set temperature
OPTIONS.tempstart = 1e1*startCost;                   % InitialTemp
OPTIONS.maxitertemp =   50*length(startGuess);       % Max-iterations per temp
OPTIONS.maxitertemp0 =  200*length(startGuess);      % Max-iterations per temp0

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

for i=1:size(OPTVAR_INDEX,2)
    dataPoints=dataPoints+length(EXPDATA.time{i});
end

CUTOFF = chi2inv(0.95,dataPoints);

file = ['parameterValues' ' ' modelName ' ' datestr(datetime('now')) '.dat'];

FID = fopen(file, 'wt');

N_iter=1;      % Number of iterations for the optimization

for i=1:N_iter
    startGuess=X;
    startCost = costFunction(startGuess);
    disp(['Initial cost with start guess: ' num2str(startCost)]);
    
    [X,FVAL,EXITFLAG] = simannealingSBAOClusteringL(@costFunction,startGuess,OPTIONS);
    
end

