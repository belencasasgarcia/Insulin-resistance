function [] = plotSimulation(param,icOrig0,plotFluxes)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global modelName
global SIMTIME
global pNames
global EXPDATA

xmin=SIMTIME(1);
xmax=SIMTIME(end);


SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig0,pNames,param);

%figure()

% Glucose in the liver and islets compartment

% Glucose in the liver compartment

subplot(2,3,1)

plot(SSsimData.time,SSsimData.variablevalues(:,3),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlim ([xmin xmax])
ylim ([0 12])

xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('Liver')

% Glucose in the islets compartment

subplot(2,3,2)

plot(SSsimData.time,SSsimData.variablevalues(:,4),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlim ([xmin xmax])
ylim ([0 12])

xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('Islets')

%Measured glucose

subplot(2,3,3)

errorbar(EXPDATA.time{1},EXPDATA.mean{1},EXPDATA.SD{1},'r.',...
    'linewidth',1)

hold on

plot(SSsimData.time,SSsimData.variablevalues(:,1),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('Measured (Liver + Islets)')

xlim ([xmin xmax])
ylim ([0 12])

% Insulin in the liver and islets compartment

% Insulin in the liver compartment

subplot(2,3,4)

plot(SSsimData.time,SSsimData.variablevalues(:,5),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('Liver')

xlim ([xmin xmax])
ylim ([0 1000])

subplot(2,3,5)

plot(SSsimData.time,SSsimData.variablevalues(:,6),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('Islets')

xlim ([xmin xmax])
ylim ([0 1000])

subplot(2,3,6)

errorbar(EXPDATA.time{2},EXPDATA.mean{2},EXPDATA.SD{2},'r.',...
    'linewidth',1)

hold on

plot(SSsimData.time,SSsimData.variablevalues(:,2),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('Measured (Liver + Islets)')

xlim ([xmin xmax])
ylim ([0 1000])

end

