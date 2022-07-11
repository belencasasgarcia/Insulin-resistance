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

load Optpara_CL

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

icOrig0=[11*param(i_V_m_hep) 11*param(i_V_m_islets) 0 0]

% Parameter changes (flow increase)
Changedparam=Optparam
%Changedparam(i_Q)=0.0002964/2;
Changedparam(i_U_ii_hep)=2*Optparam(i_U_ii_hep)


% Plot simulation with the initial parameters

figure()

plotSimulation(Optparam,icOrig0,0)

hold on

plotSimulation(Changedparam,icOrig0,0)

% Plot fluxes

figure()

plotFluxes(Optparam,icOrig0,0)

figure()

plotFluxes(Changedparam,icOrig0,0)


% Change parameters


