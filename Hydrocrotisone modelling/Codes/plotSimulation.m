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

Optparam_11mM=Optparam;
Optparam_11mM(parIndex.i_G0)=11;
Optparam_11mM(parIndex.i_I_GTT)=0;
Optparam_11mM(parIndex.i_delta_G_1)=Optparam_11mM(parIndex.i_delta_G_1_11mM);
Optparam_11mM(parIndex.i_delta_G_13)=Optparam_11mM(parIndex.i_delta_G_13_11mM);
Optparam_11mM(parIndex.i_delta_I_1)=Optparam_11mM(parIndex.i_delta_I_1_11mM);
Optparam_11mM(parIndex.i_delta_I_13)=Optparam_11mM(parIndex.i_delta_I_13_11mM);


icOrig_11mM=[(11+Optparam_11mM(parIndex.i_delta_G_1))*Optparam_11mM(parIndex.i_V_m_hep) ...
    (11+Optparam_11mM(parIndex.i_delta_G_1))*Optparam_11mM(parIndex.i_V_m_islets)...
   (0+Optparam_11mM(parIndex.i_delta_I_1))*Optparam_11mM(parIndex.i_V_m_hep) (0+Optparam_11mM(parIndex.i_delta_I_1))*Optparam(parIndex.i_V_m_islets) 0 0 ...
   5.5 0.0000000088 0 0];

SSsimData=SBPDsimulate(modelName,SIMTIME,icOrig_11mM,pNames,Optparam_11mM);

%figure()

figure()

errorbar(EXPDATA.time{1},EXPDATA.mean{1},EXPDATA.SD{1},'r.',...
    'linewidth',1)

hold on

plot(SSsimData.time,SSsimData.variablevalues(:,1),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);
xlabel ('Time (h)')
ylabel ('Glucose concentration (mM)')
title ('Liver-islet co-culture')

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
title ('Liver-islet co-culture')

% Disease progression variables

figure()

subplot(5,2,1)

plot(SSsimData.time,SSsimData.variablevalues(:,7),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Accumulated glucose (mM)')

subplot(5,2,2)

plot(SSsimData.time,SSsimData.variablevalues(:,8),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin sensitivity (^{L}/_{mU*h})')

subplot(5,2,3)

plot(SSsimData.time,SSsimData.variablevalues(:,9),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Slow glucose (mM)')

subplot(5,2,4)

plot(SSsimData.time,SSsimData.variablevalues(:,12),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Islets volume (\muL)')

subplot(5,2,5)

plot(SSsimData.time,100*SSsimData.variablevalues(:,14),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin secretion capacity (% per max)')

subplot(5,2,6)

plot(SSsimData.time,100*SSsimData.variablevalues(:,15),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Net change of beta cell volume (% per hour)')

subplot(5,2,7)

plot(SSsimData.time,100*SSsimData.variablevalues(:,29),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Reduction in insulin sensitivity due to hyperglycemia (% per max)')

subplot(5,2,8)

plot(SSsimData.time,100*SSsimData.variablevalues(:,30),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Reduction in insulin sensitivity due to an additional diabetogenic factor (% per max)')

subplot(5,2,9)

plot(SSsimData.time,100*SSsimData.variablevalues(:,29).*SSsimData.variablevalues(:,30),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Reduction in insulin sensitivity (combination) (% per max)')


%% 5.5 mM

Optparam_5p5mM=Optparam;
Optparam_5p5mM(parIndex.i_G0)=5.5;
Optparam_5p5mM(parIndex.i_I_GTT)=0;
Optparam_5p5mM(parIndex.i_delta_G_13)=Optparam_5p5mM(parIndex.i_delta_G_13_5p5mM);
Optparam_5p5mM(parIndex.i_delta_I_13)=Optparam_5p5mM(parIndex.i_delta_I_13_5p5mM);

icOrig_5p5=[(5.5)*Optparam_5p5mM(parIndex.i_V_m_hep)...
    (5.5)*Optparam_5p5mM(parIndex.i_V_m_islets) ...
    (0)*Optparam_5p5mM(parIndex.i_V_m_hep) ...
    (0)*Optparam_5p5mM(parIndex.i_V_m_islets) 0 0 5.5 0.0000000088 0 0];

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

subplot(5,2,1)

plot(SSsimData.time,SSsimData.variablevalues(:,7),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Accumulated glucose (mM)')

subplot(5,2,2)

plot(SSsimData.time,SSsimData.variablevalues(:,8),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin sensitivity (^{L}/_{mU*h})')

subplot(5,2,3)

plot(SSsimData.time,SSsimData.variablevalues(:,9),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Slow glucose (mM)')

subplot(5,2,4)

plot(SSsimData.time,SSsimData.variablevalues(:,12),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Islets volume (\muL)')

subplot(5,2,5)

plot(SSsimData.time,100*SSsimData.variablevalues(:,14),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Insulin secretion capacity (% per max)')

subplot(5,2,6)

plot(SSsimData.time,100*SSsimData.variablevalues(:,15),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Net change of beta cell volume (% per hour)')

subplot(5,2,7)

plot(SSsimData.time,100*SSsimData.variablevalues(:,29),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Reduction in insulin sensitivity due to hyperglycemia (% per max)')

subplot(5,2,8)

plot(SSsimData.time,100*SSsimData.variablevalues(:,30),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Reduction in insulin sensitivity due to an additional diabetogenic factor (% per max)')

subplot(5,2,9)

plot(SSsimData.time,100*SSsimData.variablevalues(:,29).*SSsimData.variablevalues(:,30),'Color',...
[0, 0.4470, 0.7410],'LineWidth',2);

title('Reduction in insulin sensitivity (combination) (% per max)')


end

