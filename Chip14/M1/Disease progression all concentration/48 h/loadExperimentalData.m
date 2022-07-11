%% Load and plot experimental data

%global EXPDATA_G
%global EXPDATA_I


load EXPDATA_G
load EXPDATA_I
load EXPDATA_G_hep
load EXPDATA_I_hep

%Convert input data from table to matrix 
EXPDATA_G=table2array(EXPDATA_G);
EXPDATA_I=table2array(EXPDATA_I);
EXPDATA_G_hep=table2array(EXPDATA_G_hep);
EXPDATA_I_hep=double(table2array(EXPDATA_I_hep));

% Select the first 48 hours

% Glucose data from liver+islets
EXPDATA=[];
EXPDATA.time{1} = EXPDATA_G(2:4,1);    % Time (absolute) (hours)
EXPDATA.mean{1} = EXPDATA_G(2:4,2);    % Time (absolute) (hours)
EXPDATA.SD{1} = EXPDATA_G(2:4,4);      % Time (absolute) (hours)

% Insulin data from liver+islets
EXPDATA.time{2} = EXPDATA_I(2:4,1);    % Time (absolute) (hours)
EXPDATA.mean{2} = EXPDATA_I(2:4,2);    % Time (absolute) (hours)
EXPDATA.SD{2} = EXPDATA_I(2:4,4);      % Time (absolute) (hours)

% Glucose data from liver only  
EXPDATA.time{3} = EXPDATA_G_hep(2:4,1);    % Time (absolute) (hours)
EXPDATA.mean{3} = EXPDATA_G_hep(2:4,2);    % Time (absolute) (hours)
EXPDATA.SD{3} = EXPDATA_G_hep(2:4,4);      % Time (absolute) (hours)

% Insulin data from liver only
EXPDATA.time{4} = EXPDATA_I_hep(:,1);    % Time (absolute) (hours)
EXPDATA.mean{4} = EXPDATA_I_hep(:,2);    % Time (absolute) (hours)
EXPDATA.SD{4} = EXPDATA_I_hep(:,4);      % Time (absolute) (hours)


% Insulin consumption from liver only

figure()

title ('Experimental data')

subplot(2,2,1)
errorbar(EXPDATA.time{1},EXPDATA.mean{1},EXPDATA.SD{1},'r.',...
    'linewidth',1)
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')

subplot(2,2,2)
errorbar(EXPDATA.time{2},EXPDATA.mean{2},EXPDATA.SD{2},'r.',...
    'linewidth',1)
xlabel ('Time [h]')
ylabel ('Insulin concentration (\muU/mL)')

subplot(2,2,3)
errorbar(EXPDATA.time{3},EXPDATA.mean{3},EXPDATA.SD{3},'r.',...
    'linewidth',1)
xlabel ('Time [h]')
ylabel ('Glucose concentration, spheroids (\muU/mL)')

subplot(2,2,4)
errorbar(EXPDATA.time{4},EXPDATA.mean{4},EXPDATA.SD{4},'r.',...
    'linewidth',1)
xlabel ('Time [h]')
ylabel ('Insulin concentration, spheroids (\muU/mL)')







