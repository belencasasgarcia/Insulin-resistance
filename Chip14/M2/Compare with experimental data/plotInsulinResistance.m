function [] = plotInsulinResistance(S_i,Vm,km,vh)

Int_G=[0:1:1500]    % Accumulated glucose (above 5.5 mM)

IR=(1-(Vm*Int_G.^vh)./(km^vh+Int_G.^vh));

figure()
plot(Int_G,IR)

end

