function plotSigmoid(Sigma,Alpha)

% Plot sigmoid function

V_m_islets = 0.0003;                 % (L)
V_islets = 0.0000000088;             % (L)


G=[1:1:16];

I=(V_islets*Sigma*G.^2)./(Alpha+G.^2);

figure()
plot(G,I)

xlabel ('Glucose concentration (mM)')
ylabel ('Insulin secretion rate (uU/(L*h))')

end

