# Hierarchical-Modeling-for-Porosity-Prediction

# Description
This repo implements TPCF and hierarchical modeling for porosity prediction and image reconstruction
for polar codes under AWGN channel. The implemented parts include

* Simulation
  * Model parameters estimation 
  * Model selection
  
* Case study
  * Model parameters estimation
  * Porosity prediction
  * Model selection
  * Residual diagnose
  * Regression for comparison 
  * Porosity with thresholds
  
* Reconstruction
  * Reconstruction by SA
  * Circularity

These parts are implemented by python 3 and matlab. The speed is quite acceptable for research propose. Please have fun with it !

# Requirement
os 
PIL
numpy
math
csv
matplotlib.pyplot
porespy
scipy


We recommend you to enter these directory and read the README.md in these directory to compile and install the package. 

# Usage

  * Model parameters estimation: simulationMCMH.m
  * Model selection: modelselection.m
  * Model parameters estimation: HLM_Gibbs_nopic.m
  * Porosity prediction: HLM_Gibbs_nopic.m
  * Model selection: model selection.m
  * Residual diagnose :residuals diagnose.py
  * Regression for comparison: errorplot.py
  * Porosity with thresholds: errorplot.py
  * Reconstruction by SA: reconstruction.py
  * Circularity: circularity.m


# Related paper
