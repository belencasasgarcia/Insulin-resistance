function [] = plotFluxes(param,icOrig0,plotFluxes)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global modelName
global SIMTIME
global pNames
global EXPDATA

xmin=SIMTIME(1);
xmax=SIMTIME(end);


SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig0,pNames,param);


subplot(4,3,1)

plot(SSsimData.time,SSsimData.variablevalues(:,6),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nGm,islets/Vm,islets')

subplot(4,3,2)

plot(SSsimData.time,SSsimData.variablevalues(:,7),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nGm,hep/Vm,hep')

subplot(4,3,3)

plot(SSsimData.time,SSsimData.variablevalues(:,8),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Glucose consumption')

subplot(4,3,4)

plot(SSsimData.time,SSsimData.variablevalues(:,9),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nGm,hep/Vm,hep')

subplot(4,3,5)

plot(SSsimData.time,SSsimData.variablevalues(:,10),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nGm,islets/Vm,islets')

subplot(4,3,6)

plot(SSsimData.time,SSsimData.variablevalues(:,11),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nIm,hep/Vm,islets')

subplot(4,3,7)

plot(SSsimData.time,SSsimData.variablevalues(:,12),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nIm,islets/Vm,islets')

subplot(4,3,8)

plot(SSsimData.time,SSsimData.variablevalues(:,13),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin production')

subplot(4,3,9)

plot(SSsimData.time,SSsimData.variablevalues(:,14),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*nIm,islets/Vm,hep')

subplot(4,3,10)

plot(SSsimData.time,SSsimData.variablevalues(:,15),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Q*I_m_hep/V_m_hep')

subplot(4,3,11)

plot(SSsimData.time,SSsimData.variablevalues(:,16),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin clearance')

end

