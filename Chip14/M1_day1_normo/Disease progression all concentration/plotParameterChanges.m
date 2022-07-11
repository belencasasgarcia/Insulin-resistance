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


%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData

modelName ='Islets_liver';

optModel = SBmodel(strcat(modelName,'.txt')); 
SBPDmakeMEXmodel(optModel,modelName); 
[pNames, param] = SBparameters(optModel);
icOrig0 = SBinitialconditions(optModel);

SIMTIME=[0:0.1:48];

i_V_m_hep = ismember(pNames,'V_m_hep');
i_V_hep = ismember(pNames,'V_hep');
i_V_m_islets = ismember(pNames,'V_m_islets');
i_V_islets = ismember(pNames,'V_islets');
i_Q = ismember(pNames,'Q');
i_EGP_hep = ismember(pNames,'EGP_hep');
i_U_ii_hep = ismember(pNames,'U_ii_hep');
i_S_i = ismember(pNames,'S_i');
i_Sigma = ismember(pNames,'Sigma');
i_Alpha = ismember(pNames,'Alpha');
i_CL_i_hep = ismember(pNames,'CL_i_hep');

% SBPDsimulate options
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through

file = ['parameterValues16-Sep-2019.dat'];

tmpHold = load(file);

%file = 'Optar_CL.mat'
%X=load(file)

%if isempty(tmpHold)
%    disp('Empty parameter dataset');
%else
    X = tmpHold(1,2:end);
    plotFunctionOptimal(X','Islets_liver')
%end

% Half volume in each compartment

X(i_Q)=10*X(i_Q);
plotFunctionOptimal(X,'Islets_liver')



