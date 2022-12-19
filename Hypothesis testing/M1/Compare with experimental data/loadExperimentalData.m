%% Load and plot experimental data

global EXPDATA

load EXPDATA_G_11mM_corr
load EXPDATA_I_11mM
load EXPDATA_G_hep
load EXPDATA_I_hep
load EXPDATA_G_5p5mM_corr
load EXPDATA_I_5p5mM
load EXPDATA_G_11mM_doseA
load EXPDATA_G_5p5mM_doseA
load EXPDATA_G_11mM_doseB
load EXPDATA_G_5p5mM_doseB
load EXPDATA_I_11mM_doseA
load EXPDATA_I_5p5mM_doseA
load EXPDATA_I_11mM_doseB
load EXPDATA_I_5p5mM_doseB




%Convert input data from table to matrix 
EXPDATA_G=double(table2array(EXPDATA_G_11mM));
EXPDATA_I=double(table2array(EXPDATA_I));
EXPDATA_G_hep=double(table2array(EXPDATA_G_hep));
%EXPDATA_I_hep=double(table2array(EXPDATA_I_hep));
EXPDATA_G_5p5mM=double(table2array(EXPDATA_G_5p5mM));
EXPDATA_I_5p5mM=double(table2array(EXPDATA_I_5p5mM));
EXPDATA_G_5p5mM_doseA=double(table2array(EXPDATA_G_5p5mM_doseA));
EXPDATA_G_11mM_doseA=double(table2array(EXPDATA_G_11mM_doseA));
EXPDATA_G_5p5mM_doseB=double(table2array(EXPDATA_G_5p5mM_doseB));
EXPDATA_G_11mM_doseB=double((EXPDATA_G_11mM_doseB));

EXPDATA_I_5p5mM_doseA=double(table2array(EXPDATA_I_5p5mM_doseA));
EXPDATA_I_11mM_doseA=double((EXPDATA_I_11mM_doseA));
EXPDATA_I_5p5mM_doseB=double(table2array(EXPDATA_I_5p5mM_doseB));
EXPDATA_I_11mM_doseB=double(table2array(EXPDATA_I_11mM_doseB));




% % Glucose data from liver+islets
% EXPDATA=[];
% EXPDATA.time{1} = EXPDATA_G(2:end,1);    % Time (absolute) (hours)
% EXPDATA.mean{1} = EXPDATA_G(2:end,2);    % Time (absolute) (hours)
% EXPDATA.SD{1} = EXPDATA_G(2:end,4);      % Time (absolute) (hours)
% 
% % % Insulin data from liver+islets
% EXPDATA.time{2} = EXPDATA_I(2:end,1);    % Time (absolute) (hours)
% EXPDATA.mean{2} = EXPDATA_I(2:end,2);    % Time (absolute) (hours)
% EXPDATA.SD{2} = EXPDATA_I(2:end,4);      % Time (absolute) (hours)


%% Data for modelling 

% % Glucose data from liver+islets
EXPDATA=[];
EXPDATA.time{1} = EXPDATA_G([1 2 3 4 5 6 7 8],1);    % Time (absolute) (hours)
EXPDATA.mean{1} = EXPDATA_G([1 2 3 4 5 6 7 8],2);    % Time (absolute) (hours)
EXPDATA.SD{1} = EXPDATA_G([1 2 3 4 5 6 7 8],4);      % Time (absolute) (hours)

% Insulin data from liver+islets
EXPDATA.time{2} = EXPDATA_I([1 2 3 4 5 6 7 8],1);    % Time (absolute) (hours)
EXPDATA.mean{2} = EXPDATA_I([1 2 3 4 5 6 7 8],2);    % Time (absolute) (hours)
EXPDATA.SD{2} = EXPDATA_I([1 2 3 4 5 6 7 8],4);      % Time (absolute) (hours)

% Glucose data from liver only  
EXPDATA.time{3} = EXPDATA_G_hep(3:5,1);    % Time (absolute) (hours)
EXPDATA.mean{3} = EXPDATA_G_hep(3:5,2);    % Time (absolute) (hours)
EXPDATA.SD{3} = EXPDATA_G_hep(3:5,4);      % Time (absolute) (hours)

% Insulin data from liver only
EXPDATA.time{4} = EXPDATA_I_hep(:,1);    % Time (absolute) (hours)
EXPDATA.mean{4} = EXPDATA_I_hep(:,2);    % Time (absolute) (hours)
EXPDATA.SD{4} = EXPDATA_I_hep(:,4);      % Time (absolute) (hours)

%% Glucose and insulin data 5.5 mM

% Glucose data 5.5 mM 
EXPDATA.time{5} = EXPDATA_G_5p5mM([5 6 7 8],1);    % Time (absolute) (hours)
EXPDATA.mean{5} = EXPDATA_G_5p5mM([5 6 7 8],2);    % Time (absolute) (hours)
EXPDATA.SD{5} = EXPDATA_G_5p5mM([5 6 7 8],4);      % Time (absolute) (hours)

% Insulin data 5.5 mM 
EXPDATA.time{6} = EXPDATA_I_5p5mM([5 6 7 8],1);    % Time (absolute) (hours)
EXPDATA.mean{6} = EXPDATA_I_5p5mM([5 6 7 8],2);    % Time (absolute) (hours)
EXPDATA.SD{6} = EXPDATA_I_5p5mM([5 6 7 8],4);      % Time (absolute) (hours)

% Validation data

EXPDATA.time{7} = EXPDATA_G_11mM_doseA(:,1);    % Time (absolute) (hours)
EXPDATA.mean{7} = EXPDATA_G_11mM_doseA(:,2);    % Time (absolute) (hours)
EXPDATA.SD{7} = EXPDATA_G_11mM_doseA(:,4);      % Time (absolute) (hours)

EXPDATA.time{8} = EXPDATA_G_5p5mM_doseA(:,1);    % Time (absolute) (hours)
EXPDATA.mean{8} = EXPDATA_G_5p5mM_doseA(:,2);    % Time (absolute) (hours)
EXPDATA.SD{8} = EXPDATA_G_5p5mM_doseA(:,4);      % Time (absolute) (hours)

EXPDATA.time{9} = EXPDATA_G_11mM_doseB(:,1);    % Time (absolute) (hours)
EXPDATA.mean{9} = EXPDATA_G_11mM_doseB(:,2);    % Time (absolute) (hours)
EXPDATA.SD{9} = EXPDATA_G_11mM_doseB(:,4);      % Time (absolute) (hours)

EXPDATA.time{10} = EXPDATA_G_5p5mM_doseB(:,1);    % Time (absolute) (hours)
EXPDATA.mean{10} = EXPDATA_G_5p5mM_doseB(:,2);    % Time (absolute) (hours)
EXPDATA.SD{10} = EXPDATA_G_5p5mM_doseB(:,4);      % Time (absolute) (hours)

EXPDATA.time{11} = EXPDATA_I_11mM_doseA(:,1);    % Time (absolute) (hours)
EXPDATA.mean{11} = EXPDATA_I_11mM_doseA(:,2);    % Time (absolute) (hours)
EXPDATA.SD{11} = EXPDATA_I_11mM_doseA(:,4);      % Time (absolute) (hours)

EXPDATA.time{12} = EXPDATA_I_5p5mM_doseA(:,1);    % Time (absolute) (hours)
EXPDATA.mean{12} = EXPDATA_I_5p5mM_doseA(:,2);    % Time (absolute) (hours)
EXPDATA.SD{12} = EXPDATA_I_5p5mM_doseA(:,4);      % Time (absolute) (hours)

EXPDATA.time{13} = EXPDATA_I_11mM_doseB(:,1);    % Time (absolute) (hours)
EXPDATA.mean{13} = EXPDATA_I_11mM_doseB(:,2);    % Time (absolute) (hours)
EXPDATA.SD{13} = EXPDATA_I_11mM_doseB(:,4);      % Time (absolute) (hours)

EXPDATA.time{14} = EXPDATA_I_5p5mM_doseB(:,1);    % Time (absolute) (hours)
EXPDATA.mean{14} = EXPDATA_I_5p5mM_doseB(:,2);    % Time (absolute) (hours)
EXPDATA.SD{14} = EXPDATA_I_5p5mM_doseB(:,4);      % Time (absolute) (hours)

plotExperimentalData(EXPDATA)

