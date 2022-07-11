% Main simulation file 

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

modelName ='Islets_liver_2OC';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

% figure()
% plotSimulation(param)

%% Parameter indexes

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
parIndex.i_V_sample_islets = ismember(pNames,'V_sample_islets');


% Simulation time

SIMTIME=[0:0.01:48];
time_highres=[0:0.001:48];

% Low parameter bounds

N_iter=10;      % Number of iterations for the optimization

OPTIONS.lowbounds(1)=0.09; 
OPTIONS.lowbounds(2)=1e-9;       
OPTIONS.lowbounds(3)=1e-9;
OPTIONS.lowbounds(4)=1e-9;
OPTIONS.lowbounds(5)=1e-9;

OPTIONS.lowbounds=OPTIONS.lowbounds';

% High parameter bounds

OPTIONS.highbounds(1)=10;               % U_ii_hep
OPTIONS.highbounds(2)=1e9;              % S_i
OPTIONS.highbounds(3)=1e9;              % Sigma
OPTIONS.highbounds(4)=1e9;              % Alpha
OPTIONS.highbounds(5)=80;               % CL_i_hep

OPTIONS.highbounds=OPTIONS.highbounds';

startGuess(1)=0.9;                      % U_ii_hep
startGuess(2)=0.00034;                  % S_i
startGuess(3)=2420454.55;               % Sigma
startGuess(4)=61.72;                    % Alpha
startGuess(5)=0.1;                      % CL_i_hep  

X=startGuess;
 
OPTIONS.index_Optpar=boolean(sum([i_U_ii_hep i_S_i i_Sigma ...
    i_Alpha i_CL_i_hep],2)); %Indexes of the parameters to optimize

% Fix parameters manually

Optparam=param;  %Optimal parameters

%Set parameter values manually
Optparam(i_EGP_hep)=0;
Optparam(OPTIONS.index_Optpar)=startGuess;

icOrig0=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets)...
    183*param(parIndex.i_V_m_hep) 183*param(parIndex.i_V_m_islets)];

% Plot simulation and data corresponding to the initial guess
SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig0,pNames,param);

figure()
plotSimulation(Optparam,icOrig0,0)

startCost = costFunction(startGuess);

psOpts = optimoptions(@particleswarm,'Display','iter');

% Simulated annealing
% Set temperature
OPTIONS.tempstart = 1e1*startCost;                   % InitialTemp
OPTIONS.maxitertemp =   50*length(startGuess);       % Max-iterations per temp
OPTIONS.maxitertemp0 =  200*length(startGuess);      % Max-iterations per temp0

% SBPDsimulate options
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

OPTVAR_INDEX=[1 2 3 4];    %Experimental data to optimize 

%  Cut-off value for chi-2 test

% Calculate number of data points for statistical tests
dataPoints=0;

for i=1:size(OPTVAR_INDEX,2)
    dataPoints=dataPoints+length(EXPDATA.time{i});
end

CUTOFF = chi2inv(0.95,dataPoints);

file = ['parameterValues' ' ' modelName ' ' datestr(datetime('now')) '.dat'];

FID = fopen(file, 'wt');

for i=1:N_iter
    startGuess=X';
    startCost = costFunction(startGuess);
    disp(['Initial cost with start guess: ' num2str(startCost)]);

    % Optimize
    %[optParamPS, FVALPS]=particleswarm(@costFunction, ...
    %    sum(OPTIONS.index_Optpar), OPTIONS.lowbounds,OPTIONS.lowbounds);
    
    [X,FVAL,EXITFLAG] = simannealingSBAOClusteringL(@costFunction,startGuess',OPTIONS);
    
end

% Optimal parameter vector

Optparam(OPTIONS.index_Optpar)=X(1,:);

figure()
plotSimulation(Optparam,icOrig0,0)

%Check if it matches the hepatocytes

Optparam(parIndex.i_V_islets)=0;

SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig0,pNames,Optparam);

%figure()

% Glucose when only hepatocytes are present

figure()

errorbar(EXPDATA.time{3},EXPDATA.mean{3},EXPDATA.SD{3},'r.',...
     'linewidth',1)
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('Liver only ')

hold on

plot(SSsimData.time,SSsimData.variablevalues(:,1),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

icOrig=[0.0033 0.0033 0.14652933 0.14652933];

SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig,pNames,Optparam);

%figure()

% Insulin when only hepatocytes are present 

figure()

errorbar(EXPDATA.time{4},EXPDATA.mean{4},EXPDATA.SD{4},'r.',...
     'linewidth',1)
xlabel ('Time (h)')
ylabel ('Insulin concentration (\muU/L)')
title ('Liver only ')

hold on

plot(SSsimData.time,SSsimData.variablevalues(:,2),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);








