function [] = plotSimulationMass(param,icOrig0,plotFluxes)


global modelName
global SIMTIME
global pNames
global EXPDATA

xmin=SIMTIME(1);
xmax=SIMTIME(end);

V_m_hep=param(ismember(pNames,'V_m_hep'));
V_m_islets=param(ismember(pNames,'V_m_islets'));


SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig0,pNames,param);

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
ylabel ('Insulin concentration (\muU/L)')
title ('Liver')

xlim ([xmin xmax])
ylim ([0 1000])

subplot(2,3,5)

plot(SSsimData.time,SSsimData.variablevalues(:,6),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

xlabel ('Time [h]')
ylabel ('Insulin concentration (\muU/L)')
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
ylabel ('Insulin concentration (\muU/L)')
title ('Measured (Liver + Islets)')

xlim ([xmin xmax])
ylim ([0 1000])


% Plot fluxes between different compartments

if(plotFluxes)
    
figure()

subplot(4,3,1)

plot(SSsimData.time,SSsimData.variablevalues(:,5),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nGm,islets/Vm,islets')

subplot(4,3,2)

plot(SSsimData.time,SSsimData.variablevalues(:,6),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nGm,hep/Vm,hep')

subplot(4,3,3)

plot(SSsimData.time,SSsimData.variablevalues(:,7),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Glucose consumption')

subplot(4,3,4)

plot(SSsimData.time,SSsimData.variablevalues(:,8),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nGm,hep/Vm,hep')

subplot(4,3,5)

plot(SSsimData.time,SSsimData.variablevalues(:,9),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nGm,islets/Vm,islets')

subplot(4,3,6)

plot(SSsimData.time,SSsimData.variablevalues(:,10),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nIm,hep/Vm,islets')

subplot(4,3,7)

plot(SSsimData.time,SSsimData.variablevalues(:,11),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nIm,islets/Vm,islets')

subplot(4,3,8)

plot(SSsimData.time,SSsimData.variablevalues(:,12),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin production')

subplot(4,3,9)

plot(SSsimData.time,SSsimData.variablevalues(:,13),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nIm,islets/Vm,hep')

subplot(4,3,10)

plot(SSsimData.time,SSsimData.variablevalues(:,14),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*I_m_hep/V_m_hep')

subplot(4,3,11)

plot(SSsimData.time,SSsimData.variablevalues(:,15),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin clearance')

end

