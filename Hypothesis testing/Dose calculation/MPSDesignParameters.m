% MPS design parameters

% Design parameters for the pancreas-liver 2-OC. More details on these 
% parameters can be found at:

% Casas B, Vilén L, Bauer S, Kanebratt KP, Wennberg Huldt C, et al. (2022) 
% Integrated experimental-computational analysis of a HepaRG liver-islet 
% microphysiological system for human-centric diabetes research. 
% PLOS Computational Biology 18(10): e1010587

% Scaling factor MPS - human

k_scaling_tissue=100000;

% Volumes

% Liver compartment (spheroids)

N_spheroids=40;             %Number of spheroids
V_hep_spheroid=8.50E-08;    %Volume of hepatocytes per spheroid (L)
V_hep=3.40E-06;             %Total volume of hepatocytes (L)

V_m_hep=300E-06;            %Media volume in the hepatocyte compartment (L) 

% Pancreas compartment (islets)

N_islets=10;                %Number of pancreatic islets
V_beta_islet=8.80E-10;      %Volume of beta cells per islet (L)
V_islets=8.80E-09;          %Total volume of beta cells (L)

V_m_islets=300E-06;         %Media volume in the hepatocyte compartment (L) 

%Flow rates

Q=2.96E-04;                 %Flow rate between compartment (L/h)

%Number of experimental replicates

N_replicates=3;

