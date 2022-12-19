function [] = plotSimulation(SSsimData,plotFluxes)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global modelName
global SIMTIME
global pNames
global EXPDATA
global parIndex
global Optparam
global param

xmin=SIMTIME(1);
xmax=SIMTIME(end);

%figure()

figure()

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

figure()

errorbar(EXPDATA.time{2},EXPDATA.mean{2},EXPDATA.SD{2},'r.',...
    'linewidth',1)

hold on

plot(SSsimData.time,SSsimData.variablevalues(:,2),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('Measured (Liver + Islets)')

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

% Disease progression variables

figure()

subplot(4,2,1)

plot(SSsimData.time,SSsimData.variablevalues(:,7),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Accumulated glucose (mM)')

subplot(4,2,2)

plot(SSsimData.time,SSsimData.variablevalues(:,8),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin sensitivity (^{L}/_{mU*h})')

subplot(4,2,3)

plot(SSsimData.time,SSsimData.variablevalues(:,9),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Slow glucose (mM)')

subplot(4,2,4)

plot(SSsimData.time,SSsimData.variablevalues(:,12),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Islets volume (\muL)')

subplot(4,2,5)

plot(SSsimData.time,100*SSsimData.variablevalues(:,14),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin secretion capacity (% per max)')

subplot(4,2,6)

plot(SSsimData.time,100*SSsimData.variablevalues(:,15),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Net change of beta cell volume (% per hour)')

%% 5.5 mM

Optparam_5p5mM=Optparam;
Optparam_5p5mM(parIndex.i_G0)=5.5;
Optparam_5p5mM(parIndex.i_I_GTT)=0;
Optparam_5p5mM(parIndex.i_G0_7)=5.5;
Optparam_5p5mM(parIndex.i_I0_7)=0;

icOrig_5p5=[5.5*param(parIndex.i_V_m_hep) 5.5*param(parIndex.i_V_m_islets) ...
    0*param(parIndex.i_V_m_hep) 0*param(parIndex.i_V_m_islets)...
    0 0 5.5 0.0000000088];

SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig_5p5,pNames,Optparam_5p5mM);

figure()

errorbar(EXPDATA.time{5},EXPDATA.mean{5},EXPDATA.SD{5},'r.',...
    'linewidth',1)

hold on

plot(SSsimData.time,SSsimData.variablevalues(:,1),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('5.5 mM')

xlim ([xmin xmax])
ylim ([0 12])

figure()

errorbar(EXPDATA.time{6},EXPDATA.mean{6},EXPDATA.SD{6},'r.',...
    'linewidth',1)

hold on

plot(SSsimData.time,SSsimData.variablevalues(:,2),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('5.5 mM')

figure()

subplot(4,2,1)

plot(SSsimData.time,SSsimData.variablevalues(:,7),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Accumulated glucose (mM)')

subplot(4,2,2)

plot(SSsimData.time,SSsimData.variablevalues(:,8),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin sensitivity (^{L}/_{mU*h})')

subplot(4,2,3)

plot(SSsimData.time,SSsimData.variablevalues(:,9),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Slow glucose (mM)')

subplot(4,2,4)

plot(SSsimData.time,SSsimData.variablevalues(:,12),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Islets volume (\muL)')

subplot(4,2,5)

plot(SSsimData.time,100*SSsimData.variablevalues(:,14),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin secretion capacity (% per max)')

subplot(4,2,6)

plot(SSsimData.time,100*SSsimData.variablevalues(:,15),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Net change of beta cell volume (% per hour)')
% 


end

