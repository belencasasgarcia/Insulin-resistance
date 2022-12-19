function [] = plotSimulation(SSsimData,plotFluxes,SIMTIME,modelName)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global pNames
global EXPDATA
global parIndex
global Optparam
global param

xmin=SIMTIME(1);
xmax=SIMTIME(end);

%figure()

figure()

plot(SSsimData.time,SSsimData.variablevalues(:,1),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('Measured (Liver + Islets)')

xlim ([xmin xmax])
ylim ([0 12])

figure()

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

plot(SSsimData.time,SSsimData.variablevalues(:,2),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('Measured (Liver + Islets)')

xlim ([xmin xmax])
ylim ([0 1000])


%% Plot glucose and insulin separately
%% 11 mM
% Glucose

figure()

plot(SSsimData.time,SSsimData.variablevalues(:,1),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('Measured (Liver + Islets)')

xlim ([xmin xmax])
ylim ([0 12])

figure()

plot(SSsimData.time,SSsimData.variablevalues(:,2),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('Measured (Liver + Islets)')

%% 5.5 mM

Optparam_5p5mM=Optparam;
Optparam_5p5mM(parIndex.i_G0)=5.5;
Optparam_2p8mM(parIndex.i_S_i)=0.015;

icOrig_5p5=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets) ...
    800*param(parIndex.i_V_m_hep) 800*param(parIndex.i_V_m_islets)...
    0 2e+02 5.57 0.0000000088];

SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig_5p5,pNames,Optparam_5p5mM);

figure()

plot(SSsimData.time,SSsimData.variablevalues(:,1),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('5.5 mM')

xlim ([xmin xmax])
ylim ([0 12])

figure()

plot(SSsimData.time,SSsimData.variablevalues(:,2),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('5.5 mM')

%% 2.8 mM

Optparam_2p8mM=Optparam;
Optparam_2p8mM(parIndex.i_G0)=2.8;
Optparam_2p8mM(parIndex.i_S_i)=0.015;
    
icOrig_2p8=[11*param(parIndex.i_V_m_hep) 11*param(parIndex.i_V_m_islets) ...
    800*param(parIndex.i_V_m_hep) 800*param(parIndex.i_V_m_islets)...
    0 2e+02 5.57 0.0000000088];

SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig_2p8,pNames,Optparam_2p8mM);

figure()

plot(SSsimData.time,SSsimData.variablevalues(:,1),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('2.8 mM')

xlim ([xmin xmax])
ylim ([0 12])

figure()

plot(SSsimData.time,SSsimData.variablevalues(:,2),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (mU/L)')
title ('2.8 mM')


end

