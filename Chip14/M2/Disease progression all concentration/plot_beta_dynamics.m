function [G_zero1 G_zero2] = plot_beta_dynamics(kv,d0,r1,r2)

% Plot the function describing beta cell dynamics as a function of insulin
% and calculate the zeros (no beta cell changes)
G=[0:0.01:40];

beta_dynamics=kv*(-d0+r1*G-r2*G.^2);

figure()

plot(G,beta_dynamics*100,'b')
xlabel ('Glucose (mM)')
ylabel ('Net beta cell volume change (% /hour)')

G_zero1=(r1 - (r1^2 - 4*d0*r2)^(1/2))/(2*r2);
G_zero2=(r1 + (r1^2 - 4*d0*r2)^(1/2))/(2*r2);

% Plot horizontal line corresponding to zero beta volume change
hold on
%plot(G,zeros(size(G)),'r')

% Plot vertical lines corresponding to glucose values with zero glucose
% change

line_values=[-1:0.1:5];

hold on
plot(G_zero1*ones(size(line_values)),line_values,'k')

hold on
plot(G_zero2*ones(size(line_values)),line_values,'k')

xlim([1 15])





end

