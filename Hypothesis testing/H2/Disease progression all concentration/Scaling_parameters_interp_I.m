% Change volumes

clear
clc
close all
format long
format compact

% Declare global variables

global EXPDATA_G
global EXPDATA_I
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

load Optparam_3

%% Parameter indexes

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

SIMTIME=[0:0.1:48];

% Plot simulation with the initial parameters

plotSimulation(Optparam)

% Change parameters

Optparam(i_V_m_hep)=1.65E-05/2;
Optparam(i_V_hep)=3.41E-06;
Optparam(i_V_m_islets)=1.65E-05/2;
Optparam(i_V_islets)=2.5E-07;
Optparam(i_Q)=4E-05;
%Optparam(i_Q)=1E-07;
%Optparam(i_U_ii_hep)=0;
%Optparam(i_EGP_hep)=10;

figure()
plotSimulation(Optparam)

Optparam(i_V_m_hep)=5/2;
Optparam(i_V_hep)=1.3;
Optparam(i_V_m_islets)=5/2;
Optparam(i_V_islets)=0.12;
Optparam(i_Q)=72;
