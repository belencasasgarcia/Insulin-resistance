%% Load experimental data

global EXPDATA

fd=cd;
cd('..')
cd DATA

load EXPDATA_G_11mM_corr
load EXPDATA_I_11mM
load EXPDATA_G_hep
load EXPDATA_I_hep
load EXPDATA_G_5p5mM_corr
load EXPDATA_I_5p5mM

cd(fd)

%Convert input data from table to matrix 
EXPDATA_G=double(table2array(EXPDATA_G_11mM));
EXPDATA_I=double(table2array(EXPDATA_I));
EXPDATA_G_hep=double(table2array(EXPDATA_G_hep));
EXPDATA_G_5p5mM=double(table2array(EXPDATA_G_5p5mM));
EXPDATA_I_5p5mM=double(table2array(EXPDATA_I_5p5mM));

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

% Glucose and insulin data 5.5 mM

% Glucose data 5.5 mM 
EXPDATA.time{5} = EXPDATA_G_5p5mM([5 6 7 8],1);    % Time (absolute) (hours)
EXPDATA.mean{5} = EXPDATA_G_5p5mM([5 6 7 8],2);    % Time (absolute) (hours)
EXPDATA.SD{5} = EXPDATA_G_5p5mM([5 6 7 8],4);      % Time (absolute) (hours)

% Insulin data 5.5 mM 
EXPDATA.time{6} = EXPDATA_I_5p5mM([5 6 7 8],1);    % Time (absolute) (hours)
EXPDATA.mean{6} = EXPDATA_I_5p5mM([5 6 7 8],2);    % Time (absolute) (hours)
EXPDATA.SD{6} = EXPDATA_I_5p5mM([5 6 7 8],4);      % Time (absolute) (hours)