%% Load and plot experimental data

global EXPDATA
global VALIDATIONDATA


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

% Change SD of the first dose to the original value
EXPDATA.SD{1}(1)=0.24908;
EXPDATA.SD{1}(end)=0.127954;

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

plotExperimentalData(EXPDATA)

% Data for validation

VALIDATIONDATA_G_5p5mM_doseA=double(table2array(EXPDATA_G_5p5mM_doseA));
VALIDATIONDATA_G_11mM_doseA=double(table2array(EXPDATA_G_11mM_doseA));
VALIDATIONDATA_G_5p5mM_doseB=double(table2array(EXPDATA_G_5p5mM_doseB));
VALIDATIONDATA_G_11mM_doseB=double((EXPDATA_G_11mM_doseB));

VALIDATIONDATA_I_5p5mM_doseA=double(table2array(EXPDATA_I_5p5mM_doseA));
VALIDATIONDATA_I_11mM_doseA=double((EXPDATA_I_11mM_doseA));
VALIDATIONDATA_I_5p5mM_doseB=double(table2array(EXPDATA_I_5p5mM_doseB));
VALIDATIONDATA_I_11mM_doseB=double(table2array(EXPDATA_I_11mM_doseB));


VALIDATIONDATA.time{1} = VALIDATIONDATA_G_11mM_doseA(:,1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{1} = VALIDATIONDATA_G_11mM_doseA(:,2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{1} = VALIDATIONDATA_G_11mM_doseA(:,4);      % Time (absolute) (hours)

VALIDATIONDATA.time{2} = VALIDATIONDATA_G_5p5mM_doseA(:,1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{2} = VALIDATIONDATA_G_5p5mM_doseA(:,2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{2} = VALIDATIONDATA_G_5p5mM_doseA(:,4);      % Time (absolute) (hours)

VALIDATIONDATA.time{3} = VALIDATIONDATA_G_11mM_doseB([1 3 4 5],1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{3} = VALIDATIONDATA_G_11mM_doseB([1 3 4 5],2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{3} = VALIDATIONDATA_G_11mM_doseB([1 3 4 5],5);      % Time (absolute) (hours)

VALIDATIONDATA.time{4} = VALIDATIONDATA_G_5p5mM_doseB([1 3 4 5],1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{4} = VALIDATIONDATA_G_5p5mM_doseB([1 3 4 5],2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{4} = VALIDATIONDATA_G_5p5mM_doseB([1 3 4 5],5);      % Time (absolute) (hours)

VALIDATIONDATA.time{5} = VALIDATIONDATA_I_11mM_doseA(:,1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{5} = VALIDATIONDATA_I_11mM_doseA(:,2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{5} = VALIDATIONDATA_I_11mM_doseA(:,4);      % Time (absolute) (hours)

VALIDATIONDATA.time{6} = VALIDATIONDATA_I_5p5mM_doseA(:,1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{6} = VALIDATIONDATA_I_5p5mM_doseA(:,2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{6} = VALIDATIONDATA_I_5p5mM_doseA(:,4);      % Time (absolute) (hours)

VALIDATIONDATA.time{7} = VALIDATIONDATA_I_11mM_doseB([1 3 4 5],1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{7} = VALIDATIONDATA_I_11mM_doseB([1 3 4 5],2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{7} = VALIDATIONDATA_I_11mM_doseB([1 3 4 5],4);      % Time (absolute) (hours)

VALIDATIONDATA.time{8} = VALIDATIONDATA_I_5p5mM_doseB([1 3 4 5],1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{8} = VALIDATIONDATA_I_5p5mM_doseB([1 3 4 5],2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{8} = VALIDATIONDATA_I_5p5mM_doseB([1 3 4 5],4);      % Time (absolute) (hours)

% Data corrections

% Corrections for dose B to account for offset at t=0 at both hyper and
% normoglycemia

%VALIDATIONDATA.SD{3} = VALIDATIONDATA.SD{3}+(11.017-11);     % Time (absolute) (hours)
%VALIDATIONDATA.SD{4} = VALIDATIONDATA.SD{4}+(11-9.997);      % Time (absolute) (hours)

% Corrections to account for offsets in insulin data

VALIDATIONDATA.SD{7} = VALIDATIONDATA.SD{7}+400;      % Time (absolute) (hours)
VALIDATIONDATA.SD{8} = VALIDATIONDATA.SD{8}+400;      % Time (absolute) (hours)
 
VALIDATIONDATA.SD{7}(1) = 400; 
VALIDATIONDATA.SD{8}(1) = 400;      % Time (absolute) (hours)
 
