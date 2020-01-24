# From individual behavior and pack dynamics to population responses: An individual-based approach to model the wolf social life cycle

This repository contains the files with R code to run an individual-based model (IBM) simulating wolf social life cycle including individual behavior and pack dynamics. The model is described in a paper submitted to Ecological Modelling.

The file initParam.R is to build the initial wolf population needed to launch the IBM simulations and to set the model parameters. The file submodels.R includes all sub-models used in the wolf IBM and detailed in the Methods section of the submitted paper. The file run.R runs the wolf IBM and extract some population results. It first call the submodel.R file to load all sub-models. Second, it calls the initParam.R file to create the initial population and load the model parameters. Then, a loop organizes the different sub-models, runs the simulation and records the outputs. Finally, some results are extracted from saved model outputs and figures are produced. More comments are included in each file.

Sarah Bauduin, Janvier 2020
