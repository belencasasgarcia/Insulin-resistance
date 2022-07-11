% Plot model simulations with uncertainty based on the estimated parameter
% values

clear
clc
close all
format long
format compact

plotData_H1=load('plotData_H1.mat');
plotData_H2=load('plotData_H2.mat');

plotData_H1=plotData_H1.plotData;
plotData_H2=plotData_H2.plotData;

% Define time axes to plot

simTime_5p5mM=[0:0.01:336];
simTime_11mM=[0:0.01:336];

plotFunctionUncertainty_lines(simTime_11mM, simTime_5p5mM, plotData_H1, plotData_H2)




