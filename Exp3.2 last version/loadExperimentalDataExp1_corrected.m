%% Load and plot experimental data


N_normo_noGTT=10;   %Number of biological replicates in the experiment
N=4;                %Number of biological replicates in the rest of experiments            

global EXPDATA

% Load data from MPS experiments. _corrected refers to corrections in 
% measured SD, as described in the manuscript.

fd=cd;


cd DATA

% Glucose data

load G_normo_lowHC_noGTT
load G_normo_lowHC_GTT
load G_normo_highHC_noGTT
load G_normo_highHC_GTT
load G_hyper_lowHC_GTT
load G_hyper_highHC_GTT

% Corrected glucose data
load G_hyper_highHC_GTT_corr
load G_normo_highHC_GTT_corr


% Insulin data

load I_normo_lowHC_noGTT
load I_normo_lowHC_GTT
load I_normo_highHC_noGTT
load I_normo_highHC_GTT
load I_hyper_lowHC_GTT
load I_hyper_highHC_GTT

cd(fd)

%Convert input data from table to matrix 

G_hyper_highHC_GTT=double(table2array(G_hyper_highHC_GTT));
G_hyper_lowHC_GTT=double(table2array(G_hyper_lowHC_GTT));
G_normo_highHC_GTT=double(table2array(G_normo_highHC_GTT));
I_hyper_highHC_GTT=double(table2array(I_hyper_highHC_GTT));
I_hyper_lowHC_GTT=double(table2array(I_hyper_lowHC_GTT));

G_hyper_highHC_GTT_corr=double(table2array(G_hyper_highHC_GTT_corr));
G_normo_highHC_GTT_corr=double(table2array(G_normo_highHC_GTT_corr));

%% Data for modelling 

% Estimation using first and last GTT

%% High HC


% Glucose data from hyperglycemia, high HC
EXPDATA=[];
EXPDATA.time{1} = G_hyper_highHC_GTT_corr(:,1);    % Time (absolute) (hours)
EXPDATA.mean{1} = G_hyper_highHC_GTT_corr(:,2);    % Time (absolute) (hours)
EXPDATA.SD{1} = G_hyper_highHC_GTT_corr(:,4);      % Time (absolute) (hours)

% Insulin data from hyperglycemia, high HC
EXPDATA.time{2} = I_hyper_highHC_GTT(:,1);    % Time (absolute) (hours)
EXPDATA.mean{2} = I_hyper_highHC_GTT(:,2);    % Time (absolute) (hours)
EXPDATA.SD{2} = I_hyper_highHC_GTT(:,4);      % Time (absolute) (hours)

G_normo_highHC_GTT_corr([5 6 7 8 9],4)=[2.333620898 2.089877948 2.111482674 2.093796527 2.094085896]

% Glucose data from normoglycemia, high HC, GTT
EXPDATA.time{5} = G_normo_highHC_GTT_corr([6 7 8 9],1);    % Time (absolute) (hours)
EXPDATA.mean{5} = G_normo_highHC_GTT_corr([6 7 8 9],2);    % Time (absolute) (hours)
EXPDATA.SD{5} = G_normo_highHC_GTT_corr([6 7 8 9],4);      % Time (absolute) (hours)

% Change the offset to account for the glucose level at the beginning of
% the GTT

% Change SD to account for the offset at day 13

% Insulin data from normoglycemia, high HC, GTT
EXPDATA.time{6} = I_normo_highHC_GTT([6 7],1);    % Time (absolute) (hours)
EXPDATA.mean{6} = I_normo_highHC_GTT([6 7],2);    % Time (absolute) (hours)
EXPDATA.SD{6} = I_normo_highHC_GTT([6 7],4);      % Time (absolute) (hours)

EXPDATA.SD{6}./EXPDATA.mean{6}*100


% Glucose data from normoglycemia, high HC, no GTT
% EXPDATA.time{5} = G_normo_highHC_noGTT(:,1);    % Time (absolute) (hours)
% EXPDATA.mean{5} = G_normo_highHC_noGTT(:,2);    % Time (absolute) (hours)
% EXPDATA.SD{5} = G_normo_highHC_noGTT(:,4);      % Time (absolute) (hours)
% 
% % Insulin data from normoglycemia, high HC, no GTT
% EXPDATA.time{6} = I_normo_highHC_noGTT(:,1);    % Time (absolute) (hours)
% EXPDATA.mean{6} = I_normo_highHC_noGTT(:,2);    % Time (absolute) (hours)
% EXPDATA.SD{6} = I_normo_highHC_noGTT(:,4);      % Time (absolute) (hours)

%% Low HC

% Glucose data from hyperglycemia, low HC
EXPDATA.time{7} = G_hyper_lowHC_GTT(:,1);    % Time (absolute) (hours)
EXPDATA.mean{7} = G_hyper_lowHC_GTT(:,2);    % Time (absolute) (hours)
EXPDATA.SD{7} = G_hyper_lowHC_GTT(:,4);      % Time (absolute) (hours)

% Insulin data from hyperglycemia, low HC
EXPDATA.time{8} = I_hyper_lowHC_GTT(:,1);    % Time (absolute) (hours)
EXPDATA.mean{8} = I_hyper_lowHC_GTT(:,2);    % Time (absolute) (hours)
EXPDATA.SD{8} = I_hyper_lowHC_GTT(:,4);      % Time (absolute) (hours)

% Glucose data from normoglycemia, low HC, GTT
EXPDATA.time{9} = G_normo_lowHC_GTT(:,1);    % Time (absolute) (hours)
EXPDATA.mean{9} = G_normo_lowHC_GTT(:,2);    % Time (absolute) (hours)
EXPDATA.SD{9} = G_normo_lowHC_GTT(:,4);      % Time (absolute) (hours)

% Insulin data from hyperglycemia, low HC, GTT
EXPDATA.time{10} = I_normo_lowHC_GTT(:,1);    % Time (absolute) (hours)
EXPDATA.mean{10} = I_normo_lowHC_GTT(:,2);    % Time (absolute) (hours)
EXPDATA.SD{10} = I_normo_lowHC_GTT(:,4);      % Time (absolute) (hours)

% Glucose data from normoglycemia, low HC, no GTT
% EXPDATA.time{11} = G_normo_lowHC_noGTT(:,1);    % Time (absolute) (hours)
% EXPDATA.mean{11} = G_normo_lowHC_noGTT(:,2);    % Time (absolute) (hours)
% EXPDATA.SD{11} = G_normo_lowHC_noGTT(:,4);      % Time (absolute) (hours)

% Insulin data from normoglycemia, low HC, no GTT
% EXPDATA.time{12} = I_normo_lowHC_noGTT(:,1);    % Time (absolute) (hours)
% EXPDATA.mean{12} = I_normo_lowHC_noGTT(:,2);    % Time (absolute) (hours)
% EXPDATA.SD{12} = I_normo_lowHC_noGTT(:,4);      % Time (absolute) (hours)
% % 
% 

% Correct SD in the data


EXPDATA.SD{6}=EXPDATA.mean{6}*0.1+EXPDATA.mean{6}(1);

EXPDATA.SD{9}=EXPDATA.mean{9}*0.05+EXPDATA.mean{9}(1);

EXPDATA.SD{2}(1:4)=EXPDATA.SD{2}(1:4)+EXPDATA.mean{2}(1);



% Correct SD in insulin data: The SD is about 1% of the mean value

%plotExperimentalDataExp1