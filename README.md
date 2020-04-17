# Exploring lesser-known pack dynamics mechanisms with an individual-based approach to model the wolf life cycle

This repository contains the R code to run an individual-based model (IBM) simulating wolf social life cycle including individual behavior and pack dynamics. The model is described in a paper submitted to Ecological Modelling.

The file initParam.R builds the initial wolf population needed to run the IBM simulations and set up the model parameters. The file submodels.R includes all sub-models used in the wolf IBM and detailed in the Methods section of the paper. The file run.R runs the wolf IBM and extracts some results: 1) it calls the submodel.R file to load all sub-models, 2) it calls the initParam.R file to create the initial population and load the model parameters, 3) a loop organizes the different sub-models, runs the simulation and records the outputs and 4) some results are extracted from saved model outputs and figures are produced. Detailed comments are included in each file.

Sarah Bauduin, April 2020
