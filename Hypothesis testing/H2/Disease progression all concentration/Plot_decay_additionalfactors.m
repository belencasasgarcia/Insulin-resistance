% Plot decay in insulin resistance due to hyperglycemia and additional
% diabetogenic factors

% Diabetogenic factor

% Simulation time
time_end=288;

time=[0:0.1:time_end];
%SIMTIME=[0:0.01:48];

Km_f_add=400;
Vm_f_add=1;

F_additional=(1-(Vm_f_add*time.^2)./(Km_f_add^2+time.^2));

figure()
plot(time,F_additional);