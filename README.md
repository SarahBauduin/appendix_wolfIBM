# An individual-based model to explore the impacts of lesser-known social dynamics on wolf populations

This repository contains the R code to run an individual-based model (IBM) simulating wolf social life cycle including individual behavior and pack dynamics. The model is described in a paper submitted to Ecological Modelling.

The file initParam.R builds the initial wolf population needed to run the IBM simulations and set up the model parameters. The file submodels.R includes all sub-models used in the wolf IBM and detailed in the Methods section of the paper. The file run.R runs the wolf IBM and extracts some results: 1) it calls the submodel.R file to load all sub-models, 2) it calls the initParam.R file to create the initial population and load the model parameters, 3) a loop organizes the different sub-models, runs the simulation and records the outputs and 4) some results are extracted from saved model outputs and figures are produced. Detailed comments are included in each file.

Sarah Bauduin, July 2020

Edit 29/10/2020
The package pedantics is not available on CRAN anymore. We added the last available version from the archive to install in the folder appendix_wolfIBM/R. The code line to install it is in the file submodel.R.

