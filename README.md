# Insulin-resistance

- [Repository Contents](#repository-contents)
- [Modelling and simulation](#modelling-and-simulation)
- [Visualizations](#visualizations)

## Repository Contents:

The files can be used to reproduce the results in the paper "Diseased human pancreas and liver microphysiological system for preclinical diabetes research". These files include the code to simulate the mathematical models and perform the visualizations included in the article, and the experimental data for calibration and validation of the models.

The scripts used for modelling and simulation are implemented in Matlab R2022b and use the IQM toolbox, which can be installed from: https://iqmtools.intiquan.com/

The scripts for visualization are implemented in Python 3.9.13. The datasets resulting from the simulations performed in Matlab are saved and used as inputs in the visualizations scripts 

- [**Hypothesis testing**](https://github.com/belencasasgarcia/Insulin-resistance/tree/main/Hypothesis%20testing): Folder containing Matlab codes to reproduce results in the manuscript related to hypothesis testing
- [**Hydrocortisone modelling**](https://github.com/belencasasgarcia/Insulin-resistance/tree/main/Hydrocrotisone%20modelling): Folder containing Matlab codes to reproduce results in the manuscript related to modelling of HCT effects
- [**Visualizations**](https://github.com/belencasasgarcia/Insulin-resistance/tree/main/Visualizations): Folder containing Python codes to reproduce the visualizations (plots) in the manuscript

The folders containing the Matlab codes include subfolders with the following structure:

- [**Codes**]: Matlab codes to reproduce the modelling results in the manuscript
- [**Data**]: Experimental data for modelling calibration and validation. It includes the raw data and the processed (corrected) data, as indicated by the subscript `_corr_`
- [**Models**]: Files for the mathematical models

## Modelling and simulation

### System requirements
We have tested in MATLAB R2022b on OSX Big Sur (on MacBook Air with M1 Chip and 16 GB RAM).

### Installation guide
Add all folders and subfolders of this repo to the path, by using the command `addpath(genpath('<repository-path'))`, replacing `<repository-path>` with the path to your local copy of the repo. The latest version of the IQM tools is required, which can be installed from: https://iqmtools.intiquan.com/

The typical installation time on a "normal" desktop time is ~10 min

### Demo
To reproduce the results in the manuscript, run the `plotExperimentSimulation` files within the subfolders to reproduce the corresponding results (calibration or validation) for both hypothesis testing (H1 and H2) and hydrocortisone modelling. The expected output are the Matlab figures (.fig) corresponding to Fig.2 d,e and Fig.4 a,b and Supplementary Fig.2 in the manuscript. The expected run time for generating these results on a "normal" desktop computer is ~20 min

### Instructions for use

1. Model calibration: Add the experimental data (glucose and insulin responses) to the `Data` subfolders within the `Calibration` subfolder
2. Run the `optimizationScript` script for the corresponding models (H1 and/or H2), which is included in the `Codes` subfolder. This script optimizes the parameters in the model using a simulated annealing algorithm (`simannealingSBAOClusteringL`) and the cost function defined in the `costFunction` script. The output will be the set of parameters that provide an acceptable agreement with the experimental data. The calculations for the chi-2 test are included in the `optimizationScript`
3. Run the `plotExperimentSimulation` script to visualize the model simulations (with uncertainties) against the experimental data used for calibration 
4. Run the `plotExperimentSimulation` script in the `Validation` folder to visualize the model predictions (with uncertainties) against the experimental data used for calibration.

## Visualizations

### System requirements
We have tested in Python 3.9.13 and conda 22.9.0 on OSX Big Sur (on MacBook Air with M1 Chip and 16 GB RAM).

### Installation guide
The installation guide below assume that:

1. The repository has been cloned
2. [Anaconda][conda] is installed for package and environment management

Reproduce the environment used for visualizations with the following command:

````
conda env create -f environment/conda_environment.yml
````
This creates a new environment named `manuscript_IR`. This environment can be activated by executing the following command in the terminal:

````
conda activate manuscript_IR
````

### Demo
To reproduce the results in the manuscript, run the `plotsHypothesisTesting.ipynb` and `plotsHypothesisTesting.ipynb` notebooks to generate the manuscript plots related to hypothesis testing (Figure 2) and hydrocortisone modelling (Figure 4), respectively.

[conda]: https://docs.conda.io/en/latest/
