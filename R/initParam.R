########################
# Wolf IBM, April 2020
# Model authors (in alphabetical order): Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi
#
# From the submitted paper: Exploring lesser-known pack dynamics mechanisms with an individual-based approach to model the wolf life cycle
# Paper authors: Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi, Paolo Ciucci, Christophe Duchamp, and Olivier Gimenez
# Submitted to Ecological Modelling
########################


##############################################
# Package used to build the initial population
library(NetLogoR)
##############################################


#######################################
# Initial population for our case study
# 10 packs + 5 dispersers
# 5 packs with 2 alphas (male and female, 5 years old each) and 2 pups (male and female, 1 year old each)
# 3 packs with 2 alphas (male and female, 5 years old each), 1 yearling (male, 2 years old) and 1 pup (female, 1 year old)
# 2 packs with 2 alphas (male and female, 5 years old each), 1 adult (female, 3 years old)
# 5 dispersers (3 females, 2 males of 2 years old)

# Create a data frame of the individual wolf characteristics in the initial population
initPopWolf_DF <- data.frame(
  ID = 1:43, # 43 individuals with id from 1 to 43
  sex = c(rep(c("M", "F", "M", "F"), 8), rep(c("M", "F", "F"), 2), rep("F", 3), rep("M", 2)), # "F" or "M"
  age = c(rep(c(5, 5, 1, 1), 5), rep(c(5, 5, 2, 1), 3), rep(c(5, 5, 3), 2), rep(2, 5)),
  alpha = c(rep(c(1, 1, 0, 0), 8), rep(c(1, 1, 0), 2), rep(0, 5)), # 0 = non breeder or 1 = breeder
  packID = c(rep(1:8, each = 4), rep(9:10, each = 3), rep(NA, 5)), # pack number, 10 packs of id from 1 to 10 and 5 dispersers which have NA
  disp = c(rep(0, 38), rep(1, 5)) # 0 = resident or 1 = disperser
)

# Create the initial population using the NetLogoR package
# Create an agentMatrix object
# Individual IDs are "who" in the agentMatrix and starts at 0 (automatically created when creating the individuals)
init <- function(initPopWolf_DF){
  
  # Initialize the model objects (i.e., land and wolves)
  land <- createWorld(0, 100, 0, 100) # create a fictive land (does not matter, will not be used)
  wolves <- createTurtles(n = nrow(initPopWolf_DF), world = land) # create as many wolves as nrow(initPopWolf_DF)
  wolves <- turtlesOwn(turtles = wolves, tVar = "sex", tVal = initPopWolf_DF[, "sex"]) # define each individual characteristic
  wolves <- turtlesOwn(turtles = wolves, tVar = "age", tVal = initPopWolf_DF[, "age"])
  wolves <- turtlesOwn(turtles = wolves, tVar = "alpha", tVal = initPopWolf_DF[, "alpha"])
  wolves <- turtlesOwn(turtles = wolves, tVar = "packID", tVal = initPopWolf_DF[, "packID"])
  wolves <- turtlesOwn(turtles = wolves, tVar = "disp", tVal = initPopWolf_DF[, "disp"])
  wolves <- turtlesOwn(turtles = wolves, tVar = "dismissed", tVal = 0) # 1 = dismissed from their breeding status during the current simulated year
  wolves <- turtlesOwn(turtles = wolves, tVar = "motherID", tVal = as.numeric(NA))
  wolves <- turtlesOwn(turtles = wolves, tVar = "fatherID", tVal = as.numeric(NA))
  wolves <- turtlesOwn(turtles = wolves, tVar = "cohort", tVal = as.numeric(NA))

  return(wolves)
}
#######################################


##################
# Model parameters

# Reproduction
meanPups <- 6.1 # mean litter size

# Mortality
mortalityPup <- 0.602 # mean mortality for non dispersing pups
mortalityPupSD <- 0 # sd mortality for non dispersing pups
mortalityYearling <- 0.18 # mean mortality for non dispersing yearling
mortalityYearlingSD <- 0.04 # sd mortality for non dispersing yearling
mortalityAdult <- 0.18 # mean (non density-dependent) mortality for non dispersing adult
mortalityAdultSD <- 0.04 # sd  (non density-dependent) mortality for non dispersing adult
equilibriumDens <- 30 # pack equilibrium density
terrSize <- 104 # territory size
mortalityDispPup <- 1 # mean mortality for dispersing pups
mortalityDispPupSD <- 0 # sd mortality for dispersing pups
mortalityDisp <- 0.31 # mean mortality for dispersers
mortalityDispSD <- 0 # sd mortality for dispersers

# Pack dissolvement
nAlpha1Dissolve <- 0.258 # dissolvement probability for pack with 1 breeder
nAlpha0Dissolve <- 0.846 # dissolvement probability for pack with 0 breeder
thresholdPackSize <- 4.1 # pack size threshold for dissolvement, values tested: 3.1, 4.1, 5.1

# Replacement and establishment
thresholdRelatedness <- 0.125 # relatedness threshold
probBudd <- 0.5 # probability of budding, values tested: 0.1, 0.5, 0.9

# Dispersal
meanPackSize <- 4.405 # mean pack size
sdPackSize <- 1.251 # sd pack size
probDispPup <- 0.25 # pup dispersal probability
probDispYearling <- 0.5 # yearling dispersal probability
probDispAdult <- 0.9 # adult dispersal probability

# Immigration/emigration
nImmigrants <- c(0, 1, 2) # number of immigrants entering the study area
pEmigr <- 0.1 # proportion of disperser emigrating outside of the study area

# Adoption
probAdopt <- 0.5 # probability of adoption, values tested: 0.1, 0.5, 0.9

# Running tests during the simulation to identify potential bugs
runTests <- TRUE # put FALSE for faster simulations
##################
