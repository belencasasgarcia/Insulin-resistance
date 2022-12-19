% Plot decay in spheroid functionality

t=[0:0.01:336];
km_hep=10;
Vm_hep=1;

decay=(1-(Vm_hep*t.^2)./(km_hep^2+t.^2));

figure()
plot(t,decay)