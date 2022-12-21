%% Load experimental data

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


%% Data for modelling 

% % Glucose data from liver+islets
EXPDATA=[];
EXPDATA.time{1} = EXPDATA_G([1 2 3 4 5 6 7 8],1);    % Time (absolute) (hours)
EXPDATA.mean{1} = EXPDATA_G([1 2 3 4 5 6 7 8],2);    % Time (absolute) (hours)
EXPDATA.SD{1} = EXPDATA_G([1 2 3 4 5 6 7 8],4);      % Time (absolute) (hours)

% Change SD of the first dose to the original value (prior to the data
% correction)

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

% Data for validation

% Insulin dose value
% Simulations were performed for two insulin doses: 2500 mIU/L and 3500
% mIU/L. In the publication, the second dose is included

VALIDATIONDATA_G_5p5mM_doseB=double(table2array(EXPDATA_G_5p5mM_doseB));
VALIDATIONDATA_G_11mM_doseB=double((EXPDATA_G_11mM_doseB));

VALIDATIONDATA_I_5p5mM_doseB=double(table2array(EXPDATA_I_5p5mM_doseB));
VALIDATIONDATA_I_11mM_doseB=double(table2array(EXPDATA_I_11mM_doseB));


VALIDATIONDATA.time{3} = VALIDATIONDATA_G_11mM_doseB([1 3 4 5],1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{3} = VALIDATIONDATA_G_11mM_doseB([1 3 4 5],2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{3} = VALIDATIONDATA_G_11mM_doseB([1 3 4 5],5);      % Time (absolute) (hours)

VALIDATIONDATA.time{4} = VALIDATIONDATA_G_5p5mM_doseB([1 3 4 5],1);    % Time (absolute) (hours)
VALIDATIONDATA.mean{4} = VALIDATIONDATA_G_5p5mM_doseB([1 3 4 5],2);    % Time (absolute) (hours)
VALIDATIONDATA.SD{4} = VALIDATIONDATA_G_5p5mM_doseB([1 3 4 5],5);      % Time (absolute) (hours)

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

I_0 = 400; %Insulin offset at t=0

% The insulin offset (error) at t=0 is counted 

VALIDATIONDATA.SD{7} = VALIDATIONDATA.SD{7}+I_0;      % Time (absolute) (hours)
VALIDATIONDATA.SD{8} = VALIDATIONDATA.SD{8}+I_0;      % Time (absolute) (hours)
 
VALIDATIONDATA.SD{7}(1) = I_0; 
VALIDATIONDATA.SD{8}(1) = I_0;      
