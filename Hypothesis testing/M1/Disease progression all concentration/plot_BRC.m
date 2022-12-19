% Plot parabola Topp

% d0=0.0025;           %     Death rate at zero glucose
% r1=6.3e-04;          
% r2=1.8e-05;    

d0=0.06;           %     Death rate at zero glucose
r1=0.84e-03;          
r2=0.24e-05; 

%G_fasting=[4:0.1:16];
G_fasting=[0:1:400];
rate=-d0+r1*G_fasting-r2*G_fasting.^2;

figure()
plot(G_fasting/18,rate*100)

xlabel('Glucose (mM)')
ylabel('Rate of change of beta cell mass/day (%)');
set(gca, 'fontsize', 12)
box off




%% Change time scale to hour instead

d0=0.06/24;           %     Death rate at zero glucose
r1=0.84e-03*18/24;
r1=3.565780562964245e-04;
r2=0.24e-05*18^2/24;
r2=3.475556494749250e-05;

%G_fasting=[4:0.1:16];
G_fasting=[0:1:400]/18;
rate=-d0+r1*G_fasting-r2*G_fasting.^2;

figure()
plot(G_fasting,rate*100)

%r2=0.65e-04;
G_fasting=[0:1:400]/18;
rate=-d0+r1*G_fasting-r2*G_fasting.^2;

figure()
plot(G_fasting,rate*100)

xlabel('Glucose (mM)')
ylabel('Rate of change of beta cell mass/day (%)');
set(gca, 'fontsize', 12)
box off

