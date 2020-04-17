########################
# Wolf IBM, April 2020
# Model authors (in alphabetical order): Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi
#
# From the submitted paper: Exploring lesser-known pack dynamics mechanisms with an individual-based approach to model the wolf life cycle
# Paper authors: Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi, Paolo Ciucci, Christophe Duchamp, and Olivier Gimenez
# Submitted to Ecological Modelling
########################


#######################################
# Packages used to build the sub-models
library(NetLogoR)
library(testthat)
library(pedantics)
library(SciViews)
#######################################


###########################################
# Customed functions used in the sub-models

# Sample function that works the same with 1 (length(x) == 1) or several items (length(x) >= 1) to sample
sample.vec <- function(x, ...) x[sample(length(x), ...)]

# Calculate the relatedness between each individuals in whoInd in the reference population contained in listAllInd
relatedness <- function(listAllInd, whoInd){ # allInd = list of agentMatrix, whoInd = vector of ID for wolves we want to know the degree of relatedness
  if(length(listAllInd) == 1){
    allData <- listAllInd[[1]] # combine all the outputs
  } else {
    allData <- do.call("rbind", listAllInd) # combine all the outputs
  }
  if(sum(!is.na(allData@.Data[,"fatherID"])) == 0 & sum(!is.na(allData@.Data[,"fatherID"])) == 0){ # no info for all individuals on their mother and father
    matrixNA <- matrix(nrow = length(whoInd), ncol = length(whoInd), data = 0) # give 0 of relatedness (no related) if no info
    colnames(matrixNA) <- whoInd
    rownames(matrixNA) <- whoInd
    return(matrixNA)
  } else {
    allDataPed <- as.data.frame(cbind(allData@.Data[,"who"], allData@.Data[,"motherID"], allData@.Data[,"fatherID"], allData@.Data[,"cohort"]))
    names(allDataPed) <- c("id", "dam", "sire", "cohort") # dam = mother/female, sire = father/male, names necessary for pedantics
    allDataPedUnique <- unique(allDataPed) # remove duplicates = individuals alive more than a year (combination of ID, motherID, fatherID and cohort more than once)
    if(length(unique(allDataPedUnique[, "cohort"][!is.na(allDataPedUnique[, "cohort"])])) == 1){ # if there is only one cohort
      wolvesPedigreeSummary <- pedigreeStats(allDataPedUnique[, 1:3], graphicalReport = "n") # do not include it
    } else {
      wolvesPedigreeSummary <- pedigreeStats(allDataPedUnique[, 1:3], cohorts = allDataPedUnique[, "cohort"], graphicalReport = "n")
    }
    relatednessMatrix <- wolvesPedigreeSummary$Amatrix # matrix of relatedness coefficient between all pairs of individuals
    return(relatednessMatrix[as.numeric(rownames(wolvesPedigreeSummary$Amatrix)) %in% whoInd, 
                             as.numeric(colnames(wolvesPedigreeSummary$Amatrix)) %in% whoInd])
  }
}

# Density dependent mortality from Cubaynes et al. 2014 (Fig. 3)
mortalityDD <- function(popDens){ # popDens in km2
  popDens1000 <- popDens * 1000 # needs to be per 1000 km2
  standPop <- (popDens1000 - 53.833) / 17.984 # standardize the population density with mean and sd from Cubaynes et al. 2014
  logitPhi <- 1.196 + (-0.505 * standPop) # intercept and slope from Cubaynes et al. 2014
  phi <- 1 / (1 + exp(-logitPhi)) # back transform the logit into the mortality probability
  pMort <- 1 - phi # phi is survival, we need mortality
  return(pMort)
}
###########################################


############
# Submodels
############

##############
# Reproduction
repro <- function(wolves, allWolvesID, yearSim, popSim){
  
  # Identify the breeding individuals
  alphaInd <- NLwith(agents = wolves, var = "alpha", val = 1)
  # Breeding females
  alphaFemale <- NLwith(agents = alphaInd, var = "sex", val = "F")
  packAlphaFemale <- of(agents = alphaFemale, var = "packID")
  if(runTests){
    expect_true(all(table(packAlphaFemale) <= 1)) # only one breeding female per pack
  }
  # Breeding males
  alphaMale <- NLwith(agents = alphaInd, var = "sex", val = "M")
  packAlphaMale <- of(agents = alphaMale, var = "packID")
  if(runTests){
    expect_true(all(table(packAlphaMale) <= 1)) # only one breeding male per pack
  }
  
  # packID where there is a breeding pair
  packReproduce <- intersect(packAlphaFemale, packAlphaMale)
  # Females that will reproduce
  femaleReproduce <- NLwith(agents = alphaFemale, var = "packID", val = packReproduce)
  
  if(NLcount(femaleReproduce) != 0){
    
    if(runTests){
      numWolves <- NLcount(wolves)
    }
    
    IDFemaleReproduce <- of(agents = femaleReproduce, var = "who")
    # Number of pups each female will have (different number per female)
    nPups <- rpois(n = length(IDFemaleReproduce), lambda = meanPups)
    
    if(sum(nPups) != 0){
      
      # Remove from the loop the females which have 0 pup
      IDFemaleReproduce <- IDFemaleReproduce[nPups != 0] 
      nPups <- nPups[nPups != 0]
      femaleReproduce <- turtle(turtles = wolves, who = IDFemaleReproduce)
      
      # Create the new wolves with hatch()
      wolves <- hatch(turtles = wolves, who = IDFemaleReproduce, n = nPups, breed = "newborn") # breed = "newborn" to recognize the pups in the wolves object
      newborn <- NLwith(agents = wolves, var = "breed", val = "newborn")
      # The newborns inherit all the parent's (femaleReproduce) parameters except for the who numbers
      # Some of the inherited variables must be changed
      # The who numbers (IDs) also need to be updated so that newborns never have an IDs of an already dead wolf from this population (this is needed to later define the pedigree)
      uniqueWho <- seq(from = max(allWolvesID) + 1, to = max(allWolvesID) + NLcount(newborn), by = 1)
      
      # IDs of the fathers need to be ordered to match the corresponding mother
      maleReproduce <- NLwith(agents = alphaMale, var = "packID", val = packReproduce)
      maleReproduceData <- of(agents = maleReproduce, var = c("who", "packID"))
      IDMaleReproduce <- maleReproduceData[match(of(agents = femaleReproduce, var = "packID"), maleReproduceData[, "packID"]), "who"]
      
      # Update the new wolves variables
      wolves <- NLset(turtles = wolves, agents = newborn, 
                      var = c("sex", "age", "alpha", "breed", "who", "motherID", "fatherID", "cohort"),
                      val = data.frame(sex = sample(c("F", "M"), NLcount(newborn), replace = TRUE),
                                       age = 0,
                                       alpha = 0,
                                       breed = "turtle",
                                       who = uniqueWho,
                                       motherID = rep(IDFemaleReproduce, nPups),
                                       fatherID = rep(IDMaleReproduce, nPups),
                                       cohort = yearSim))
      
      # Update allWolvesID with the new wolves ID
      allWolvesID <- c(allWolvesID, uniqueWho)
      
      if(runTests){
        # Make sure pups were integrated in the wolves object
        expect_equivalent(NLcount(wolves), numWolves + sum(nPups)) 
      }
    }
  }
  
  if(runTests){
    # There should not be duplicated IDs
    expect_equal(length(allWolvesID), length(unique(allWolvesID)))
  }
  
  # After reproduction and new individuals are in the population, calculate the relatedness between all wolves
  temporaryOutputs <- popSim
  temporaryOutputs[[length(popSim) + 1]] <- wolves
  allWolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = of(agents = wolves, var = "who"))
  
  return(list(wolves, allWolvesID, allWolvesRelatedness))
}
##############

#######
# Aging
aging <- function(wolves){
  
  ageWolf <- of(agents = wolves, var = "age")
  
  if(runTests){
    # No wolf should be older than 15 years old
    expect_true(all(ageWolf <= 15))
  }
  
  # All wolves get 1 year older
  wolves <- NLset(turtles = wolves, agents = wolves, var = "age", val = ageWolf + 1) 
  
  return(wolves)
}
#######

###########
# Mortality
mortality <- function(wolves){
  
  if(runTests){
    numWolves <- NLcount(wolves) 
  }
  
  # Calculate the current population density before any mortality event
  # The density does not include the pups of the year (Cubaynes et al. 2014)
  withoutPups  <- NLwith(agents = wolves, var = "age", val = 2:16)
  popDens <- NLcount(withoutPups) / (CarryingCapacity * terrSize)
  # Number of packs in the population
  numPacks <- unique(of(agents = wolves, var = "packID"))
  numPacks <- numPacks[!is.na(numPacks)]
  
  # Kill the too old wolves
  wolvesOld <- NLwith(agents = wolves, var = "age", val = 16)
  wolves <- die(turtles = wolves, who = of(agents = wolvesOld, var = "who"))
  
  # Pup mortality
  wolvesPup <- NLwith(agents = wolves, var = "age", val = 1)
  if(runTests){
    # There should not be any dispersing pup at that moment
    expect_true(NLcount(NLwith(agents = wolvesPup, var = "disp", val = 1)) == 0)
  }
  IDPup <- of(agents = wolvesPup, var = "who")
  deadPup <- rbinom(n = length(IDPup), size = 1,
                    prob = rnorm(1, mean = mortalityPup, sd = mortalityPupSD)) # 0 = survive, 1 = die
  wolves <- die(turtles = wolves, who = IDPup[deadPup == 1])
  
  # Yearling mortality
  # Non-dispersing yearlings
  wolvesYearling <- NLwith(agents = wolves, var = "age", val = 2)
  nonDispYearling <- NLwith(agents = wolvesYearling, var = "disp", val = 0)
  IDnonDispYearling <- of(agents = nonDispYearling, var = "who")
  deadYearling <- rbinom(n = length(IDnonDispYearling), size = 1,
                         prob = rnorm(1, mean = mortalityYearling, sd = mortalityYearlingSD)) # 0 = survive, 1 = die
  wolves <- die(turtles = wolves, who = IDnonDispYearling[deadYearling == 1])
  # Yearlings that dispersed as pups and did not become adoptee (i.e., still dispersers) 
  Disperser <- NLwith(agents = wolves, var = "disp", val = 1)
  DisperserYearling <- NLwith(agents = Disperser, var = "age", val = 2)
  IDDisperserYearling <- of(agents = DisperserYearling, var = "who")
  deadDisperserYearling <- rbinom(n = length(IDDisperserYearling), size = 1,
                                  prob = rnorm(1, mean = mortalityDispPup, sd = mortalityDispPupSD)) # 0 = survive, 1 = die
  wolves <- die(turtles = wolves, who = IDDisperserYearling[deadDisperserYearling == 1])
  
  # Adult mortality
  # Non-dispersing adults
  wolvesAdult <- NLwith(agents = wolves, var = "age", val = 3:15)
  nonDispAdult <- NLwith(agents = wolvesAdult, var = "disp", val = 0)
  IDnonDispAdult <- of(agents = nonDispAdult, var = "who")
  if(length(numPacks) == CarryingCapacity){ # at carrying capacity, mortality density dependent
    deadAdult <- rbinom(n = length(IDnonDispAdult), size = 1,
                        prob = rnorm(1, mean = mortalityDD(popDens = popDens), sd = 0))
  } else { # there are less packs than the maximum allowed for this study area, mortality not density dependent
    deadAdult <- rbinom(n = length(IDnonDispAdult), size = 1,
                        prob = rnorm(1, mean = mortalityAdult, sd = mortalityAdultSD))
  }
  wolves <- die(turtles = wolves, who = IDnonDispAdult[deadAdult == 1])
  # Dispersersing adults
  DisperserAdult <- NLwith(agents = Disperser, var = "age", val = 3:15)
  IDDisperserAdult <- of(agents = DisperserAdult, var = "who")
  deadDisperserAdult <- rbinom(n = length(IDDisperserAdult), size = 1,
                               prob = rnorm(1, mean = mortalityDisp, sd = mortalityDispSD)) # 0 = survive, 1 = die
  wolves <- die(turtles = wolves, who = IDDisperserAdult[deadDisperserAdult == 1])
  
  if(runTests){
    # All dead individuals were removed from wolves
    expect_equal(NLcount(wolves),
                 numWolves - NLcount(wolvesOld) - sum(deadPup) - sum(deadYearling) - sum(deadAdult) - sum(deadDisperserYearling) - sum(deadDisperserAdult))
  }
  
  return(wolves)
}
###########

###################
# Pack dissolution
PackDissolvement <- function(wolves){
  
  # Dissolution of the packs where there are only 1-year old wolves
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesMature <- NLwith(agents = wolves, var = "age", val = 2:15) # packID where there is at least one non-pup wolf
  packIDwolvesMature <- unique(of(agents = wolvesMature, var = "packID"))
  PackOnlyWithPup <- wolvesPackID[!wolvesPackID %in% packIDwolvesMature] # packID with only pups in it
  PackOnlyWithPup <- PackOnlyWithPup[!is.na(PackOnlyWithPup)] # remove the NA from the dispersers
  pupDissolve <- NLwith(agents = wolves, var = "packID", val = PackOnlyWithPup)
  # Dissolve the concerned packs
  wolves <- NLset(turtles = wolves, agents = pupDissolve,
                  var = c("packID", "disp"), val = cbind(rep(NA, NLcount(pupDissolve)), rep(1, NLcount(pupDissolve))))

  if(runTests){
    # There must not be packID related only to 1 year old wolves
    newPackID <- unique(of(agents = wolves, var = "packID"))
    wolvesMature <- NLwith(agents = wolves, var = "age", val = 2:15)
    packIDwolvesMature <- unique(of(agents = wolvesMature, var = "packID"))
    expect_true(setequal(newPackID[!is.na(newPackID)], packIDwolvesMature[!is.na(packIDwolvesMature)])) #same elements but can be in different orders
  }
  
  # Possible pack dissolution regarding how many breeding individuals there are left
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  totalPackDiss <- 0 # record the number of pack that dissolve
  
  for(eachPack in wolvesPackID){
    
    # Calculate pack size
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    packSize <- NLcount(wolvesInPack)
    
    # Only small pack can dissolve
    if(packSize < thresholdPackSize){ 
      
      # How many breeding individuals there are in the pack
      nAlpha <- NLcount(NLwith(agents = wolvesInPack, var = "alpha", val = 1))
      if(nAlpha != 2){
        if(nAlpha == 1){ # only 1 breeder leaft
          packDissolve <- rbinom(n = 1, size = 1, prob = nAlpha1Dissolve)
        } else if(nAlpha == 0){ # no breeder left
          packDissolve <- rbinom(n = 1, size = 1, prob = nAlpha0Dissolve)
        }
        if(packDissolve == 1){
          totalPackDiss <- totalPackDiss + 1
          # Dissolve the concerned packs
          wolves <- NLset(turtles = wolves, agents = wolvesInPack,
                          var = c("alpha", "packID", "disp"),
                          val = cbind(alpha = rep(0, packSize), rep(NA, packSize), rep(1, packSize)))
        }
      }
    }
  }
  
  if(runTests){
    # All dissolved packs should be removed from wolves
    newWolvesPackID <- unique(of(agents = wolves, var = "packID"))
    newWolvesPackID <- newWolvesPackID[!is.na(newWolvesPackID)] # remove the NA from the dispersers
    expect_equal(length(newWolvesPackID), length(wolvesPackID) - totalPackDiss)
  }
  
  return(wolves)
}
###################

#################################################
# Replacement of breeding females by subordinates
FemaleAlphaSurbordinateReplacement <- function(wolves, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  for(eachPack in wolvesPackID){
    
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack) # wolves in the pack
    FemaleInPack <- NLwith(agents = wolvesInPack, var = "sex", val = "F") # females in the pack
    FemaleAlphaInPack <- NLwith(agents = FemaleInPack, var = "alpha", val = 1) # breeding female in the pack
    
    # If there is no breeding female but other females available in the pack
    if(NLcount(FemaleAlphaInPack) == 0 & NLcount(FemaleInPack) != 0){
      FemaleSMatureAvailable <- NLwith(agents = FemaleInPack, var = "age", val = c(2:15)) # select among these subordinates the mature ones
      IDFemaleSMatureAvailable <- as.numeric(of(agents = FemaleSMatureAvailable, var = "who")) # select their ID
      IDNewAlphaF <- sample.vec(IDFemaleSMatureAvailable, 1, replace = FALSE) # select randomly one female among these selected females
      # Change the breeding status of the new breeding female
      wolves <- NLset(turtles = wolves,
                      agents = turtle(turtles = wolves, who = IDNewAlphaF),
                      var = "alpha", val = 1)
      
      # Check if the breeding male of the pack is not related to the new breeding female
      maleAlphaInPack <- NLwith(agents = wolvesInPack, var = "alpha", val = 1) # identify the breeding male (i.e., the only breeder as there is no breeding female)
      if(length(IDNewAlphaF) != 0 & NLcount(maleAlphaInPack) != 0){
        IDmaleAlphaInPack <- of(agents = maleAlphaInPack, var = "who")
        # Relatedness between the breeders
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDNewAlphaF,
                                                 colnames(allWolvesRelatedness) == IDmaleAlphaInPack]
        # If the breeders are too related, the male looses his breeding position and becomes subordinate
        if(relatAlphaCouple > thresholdRelatedness){ 
          wolves <- NLset(turtles = wolves, agents = maleAlphaInPack, var = c("alpha", "dismissed"),
                          val = cbind(alpha = 0, dismissed = 1))
        }
      }
    }
  }
  
  if(runTests){
    # All packs with mature females (>= 2 yrs old) must have a breeding female now
    # So there must be as many unique packID with mature females in it, than the number of alpha females in these packs
    allFem <- NLwith(agents = wolves, var = "sex", val = "F")
    allFemInPack <- NLwith(agents = allFem, var = "disp", val = 0)
    matureFemInPack <- NLwith(agents = allFemInPack, var = "age", val = 2:15)
    expect_equal(length(unique(of(agents = matureFemInPack, var = "packID"))),
                 NLcount(NLwith(agents = matureFemInPack, var = "alpha", val = 1)))
  }
  
  return(wolves)
}
##################################################

###########
# Dispersal
dispersal <- function(wolves){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  # Keep in memory the original number of dispersers
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) 
  totalDispCreated <- 0
  
  # Prepare the vector to store the maximum size for each pack
  max_sizes <- numeric(0)
  if(length(wolvesPackID) != 0){ 
    max_sizes <- rep(NA, max(wolvesPackID)) # to be returned at the end of the function to be used in another function
  }
  
  #Identify the individuals that need to disperse
  for(eachPack in wolvesPackID){
    
    # How many individuals need to disperse in this pack
    packSize <- NLcount(NLwith(agents = wolves, var = "packID", val = eachPack))
    # Define the maximum pack size for this pack
    maxPackSize <- round(rnorm(n = 1, mean = meanPackSize, sd = sdPackSize))
    max_sizes[eachPack] <- maxPackSize # fill the vector
    numNeedDisp <- packSize - maxPackSize
    
    # If the pack creates dispersers
    if(numNeedDisp > 0){
      # Identify the dispersers
      wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack) # wolves in the pack
      potentialDisp <- NLwith(agents = wolvesInPack, var = "alpha", val = 0) # breeders cannot disperse
      potentialDispPup <- of(agents = NLwith(agents = potentialDisp, var = "age", val = 1), var = "who")
      potentialDispYearling <- of(agents = NLwith(agents = potentialDisp, var = "age", val = 2), var = "who")
      potentialDispAdult <- of(agents = NLwith(agents = potentialDisp, var = "age", val = 3:15), var = "who")
      IDwhichDisp <- c(potentialDispPup, potentialDispYearling, potentialDispAdult)
      
      # If in the pack there are more potential dispersers than individuals that need to disperse, select them based on their dispersal probability related to their age
      if(length(IDwhichDisp) > numNeedDisp){
        sumProb <- (length(potentialDispPup) * probDispPup) + (length(potentialDispYearling) * probDispYearling) + (length(potentialDispAdult) * probDispAdult)
        probPotentialDispPup <- rep(probDispPup / sumProb, length(potentialDispPup))
        probPotentialDispYearling <- rep(probDispYearling  / sumProb, length(potentialDispYearling))
        probPotentialDispAdult <- rep(probDispAdult  / sumProb, length(potentialDispAdult))
        # Select as many dispersers as needed according to their probabilities
        IDwhichDisp <- sample(IDwhichDisp, size = numNeedDisp, replace = FALSE, prob = c(probPotentialDispPup, probPotentialDispYearling, probPotentialDispAdult))
      }
      
      # Update the selected dispersers
      totalDispCreated <- totalDispCreated + length(IDwhichDisp) # increment the number of dispersers created
      wolves <- NLset(turtles = wolves,
                      agents = turtle(turtles = wolves, who = IDwhichDisp),
                      var = c("packID", "disp"), val = cbind(packID = rep(NA, length(IDwhichDisp)), disp = rep(1, length(IDwhichDisp))))
    }
  }
  
  if(runTests){
    # Check if all dispersers have been updated in wolves
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) # new number of dispersers
    expect_equal(numDisp + totalDispCreated, newNumDisp)
  }
  
  return(list(wolves, max_sizes))
}
###########

#############
# Immigration
immigration <- function(wolves, allWolvesID, popSim){
  
  # How many migrants are coming
  nMigrInd <- sample.vec(nImmigrants, 1, replace = FALSE)
  
  if(nMigrInd != 0){
    
    numWolves <- NLcount(wolves) # number of wolves in the population before the arrival of the immigrats
    
    # Create migrants with the same characteristicts as the wolf
    migrants <- createTurtles(n = nMigrInd, coords = wolves[1]@.Data[, c("xcor", "ycor"), drop = FALSE]) # locate them where was the first wolf
    migrants <- turtlesOwn(turtles = migrants, tVar = "sex", tVal = sample.vec(c("F", "M"), nMigrInd, replace = TRUE))
    ageMigrants <- rpois(n = nMigrInd, lambda = 2) # migrants are more likely to be young dispersers (yearlings)
    ageMigrants[ageMigrants < 1] <- 1 # pups have their age = 1
    ageMigrants[ageMigrants > 15] <- 15 # wolf age limit
    migrants <- turtlesOwn(turtles = migrants, tVar = "age", tVal = ageMigrants)
    migrants <- turtlesOwn(turtles = migrants, tVar = "alpha", tVal = 0)
    migrants <- turtlesOwn(turtles = migrants, tVar = "packID", tVal = as.numeric(NA))
    migrants <- turtlesOwn(turtles = migrants, tVar = "disp", tVal = 1)
    migrants <- turtlesOwn(turtles = migrants, tVar = "dismissed", tVal = 0)
    migrants <- turtlesOwn(turtles = migrants, tVar = "motherID", tVal = as.numeric(NA))
    migrants <- turtlesOwn(turtles = migrants, tVar = "fatherID", tVal = as.numeric(NA))
    migrants <- turtlesOwn(turtles = migrants, tVar = "cohort", tVal = as.numeric(NA))

    # Change the who number (ID) of the migrants to be unique among all the wolf ever created in this population
    uniqueWho <- seq(from = max(allWolvesID) + 1, to = max(allWolvesID) + nMigrInd, by = 1)
    migrants <- NLset(turtles = migrants, agents = migrants, var = "who", val = uniqueWho)
    # Update allWolvesID with the new wolves
    allWolvesID <- c(allWolvesID, uniqueWho)
    
    # Join the migrants to the wolf population
    wolves <- turtleSet(wolves, migrants)

    if(runTests){
      expect_equivalent(NLcount(wolves), numWolves + nMigrInd) # to make sure migrants were integrated inside the wolves object
      expect_equal(length(allWolvesID), length(unique(allWolvesID))) # there should not be duplicated IDs
    }
    
    # After immigration, new wolves have arrived, calculate the relatedness of all wolf pairs
    temporaryOutputs <- popSim
    temporaryOutputs[[length(popSim) + 1]] <- wolves
    allWolvesRelatedness <- relatedness(listAllInd = temporaryOutputs, whoInd = of(agents = wolves, var = "who"))
    
  }
  
  return(list(wolves, allWolvesID, allWolvesRelatedness))
}
#############

############
# Emigration
emigration <- function(wolves){
  
  # How many dispersers will emigrate
  dispInd <- NLwith(agents = wolves, var = "disp", val = 1) # dispersers
  nEmigr <- round(pEmigr * NLcount(dispInd))
  
  if(nEmigr != 0){
    emigrInd <- nOf(agents = dispInd, n = nEmigr) # select the emigrants
    # Kill the emigrants to remove them from the population
    wolves <- die(turtles = wolves, who = of(agents = emigrInd, var = "who"))
  }
  
  return(wolves)
}
############

##########
# Adoption
Adoptee <-  function(wolves, max_sizes){ 
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  # Shuffle wolvesPackID so that it's not always the packs at the beginning (small packID) that receveive adoptees
  wolvesPackID <- sample.vec(wolvesPackID)
  
  # Keep the original number of dispersers to test for later
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
  totalAdoptee <- 0
  
  for(eachPack in wolvesPackID){
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    packSize <- NLcount(wolvesInPack)
    
    # The packs needs to be able to adopt (i.e., have available places and regarding the probability)
    if(packSize < max_sizes[eachPack] & rbinom(n = 1, size = 1, prob = probAdopt)){ 
      
      PlacesAvailable <- max_sizes[eachPack] - packSize # how many adoptees the pack can adopt
      Disperser <- NLwith(agents = wolves, var = "disp", val = 1) # current dispersers
      PotentialAdoptee <- NLwith(agents = Disperser, var = "age", val = c(1, 2, 3)) # dispersers of 1, 2 or 3 years old
      PotentialAdopteeM <- NLwith(agents = PotentialAdoptee, var = "sex", val = "M") # male dipersers are favored for adoption
      PotentialAdopteeF <- NLwith(agents = PotentialAdoptee, var = "sex", val = "F")
      IDPotentialAdopteeM  <- of(agents = PotentialAdopteeM, var = "who") #select their ID
      IDPotentialAdopteeF  <- of(agents = PotentialAdopteeF, var = "who") #select their ID
      
      # The probability of adopting is only driven by probAdopt but not density-dependent
      # First select the adoptees among the male dispersers
      if(length(IDPotentialAdopteeM) > PlacesAvailable){
        IDRealAdoptee <- sample.vec(IDPotentialAdopteeM, PlacesAvailable, replace = FALSE)  # pick potential adoptees among the available places
      } else {
        IDRealAdoptee <- IDPotentialAdopteeM
      }
      # After choosing the males which are favored, check if there are still places in the pack to adopt females
      if(length(IDRealAdoptee) < PlacesAvailable){
        if(length(IDPotentialAdopteeF) > (PlacesAvailable - length(IDRealAdoptee))){
          IDRealAdoptee <- c(IDRealAdoptee, sample.vec(IDPotentialAdopteeF, (PlacesAvailable - length(IDRealAdoptee)), replace = FALSE))  # pick potential adoptees among the available places
        } else {
          IDRealAdoptee <- c(IDRealAdoptee, IDPotentialAdopteeF)
        }
      }
      
      # Update the adoptees
      totalAdoptee <- totalAdoptee + length(IDRealAdoptee)
      wolves <- NLset(turtles = wolves,
                      agents = turtle(turtles = wolves, who = IDRealAdoptee),
                      var = c("packID", "disp"), val = cbind(packID = rep(eachPack, length(IDRealAdoptee)),
                                                             disp = rep(0, length(IDRealAdoptee))))
    }
  }
  
  if(runTests){
    # The former dispersers that were adopted are no longer dispersers
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) # new number of dispersers
    expect_equal(numDisp - totalAdoptee, newNumDisp)
  }
  
  return(wolves)
}
##########

#######################################
# Replacement of breeders by dispersers
AlphaDisperserReplacement <- function(wolves, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  # Shuffle wolvesPackID so that it's not always the packs at the beginning (small packID) that receveive dispersers
  wolvesPackID <- sample.vec(wolvesPackID)
  
  if(runTests){
    numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1)) # keep the original number of dispersers to test for later
    numAlpha <- NLcount(NLwith(agents = wolves, var = "alpha", val = 1)) # and the number of breeders
  }
  totalReplace <- 0
  
  for(eachPack in wolvesPackID){
    
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    AlphaInPack <- NLwith(agents = wolvesInPack, var = "alpha", val = 1) # breeders in the pack
    Disperser <- NLwith(agents = wolves, var = "disp", val = 1) # all dispersers in the population

    # Breeding female replacement
    FemaleAlphaInPack <- NLwith(agents = AlphaInPack, var = "sex", val = "F") # breeding female in the pack
    if(NLcount(FemaleAlphaInPack) == 0){ # if there is no breeding female in the pack
      FemaleDisperserAvailable <- NLwith(agents = Disperser, var = "sex", val = "F") # select the female dispersers
      FemaleDMatureAvailable <- NLwith(agents = FemaleDisperserAvailable, var = "age", val = c(2:15)) # select the mature ones
      IDFemaleDMatureAvailable <- of(agents = FemaleDMatureAvailable, var = "who") # select their ID
      
      # If there is a breeding male in this pack, remove from IDFemaleDMatureAvailable the females that are too closely related
      MaleAlphaInPack <- NLwith(agents = AlphaInPack, var = "sex", val = "M") # breeding male in the pack
      if(NLcount(MaleAlphaInPack) == 1 & NLcount(FemaleDMatureAvailable) != 0){
        IDmaleAlphaInPack <- of(agents = MaleAlphaInPack, var = "who")
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDmaleAlphaInPack,
                                                 colnames(allWolvesRelatedness) %in% IDFemaleDMatureAvailable, drop = FALSE]
        relatedFemales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of females too closely related
        IDFemaleDMatureAvailable <- IDFemaleDMatureAvailable[!IDFemaleDMatureAvailable %in% relatedFemales] # remove them from the potential females to become breeder
      }
      
      # Select randomly one female among the selected females
      IDNewAlphaF <- sample.vec(IDFemaleDMatureAvailable, 1, replace = FALSE) 
      # Replace the missing breeding female by the selected one
      wolves <- NLset(turtles = wolves,
                      agents = turtle(turtles = wolves, who = IDNewAlphaF),
                      var = c("alpha", "packID", "disp"), val = cbind(alpha = 1, packID = eachPack, disp = 0))
      
      if(length(IDNewAlphaF) == 1){ # update the counter
        totalReplace <- totalReplace + 1
      }
    }
    
    # Breeding male replacement
    MaleAlphaInPack <- NLwith(agents = AlphaInPack, var = "sex", val = "M")
    if(NLcount(MaleAlphaInPack) == 0){ # if there is no breeding male in the pack
      MaleDisperserAvailable <- NLwith(agents = Disperser, var = "sex", val = "M") # select the male dispersers
      MaleDMatureAvailable <- NLwith(agents = MaleDisperserAvailable, var = "age", val = c(2:15)) # select the mature ones
      IDMaleDMatureAvailable <- of(agents = MaleDMatureAvailable, var = "who") # select their ID
      
      # If there is a breeding female in this pack, remove from IDMaleDMatureAvailable the males that are too closely related
      FemaleAlphaInPack <- NLwith(agents = AlphaInPack, var = "sex", val = "F") # breeding female already in the pack (Warning: does not include the disperser that just became breeder)
      if(NLcount(FemaleAlphaInPack) == 1 & NLcount(MaleDMatureAvailable) != 0){
        IDfemaleAlphaInPack <- of(agents = FemaleAlphaInPack, var = "who")
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDfemaleAlphaInPack,
                                                 colnames(allWolvesRelatedness) %in% IDMaleDMatureAvailable, drop = FALSE]
        relatedMales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of male too closely related
        IDMaleDMatureAvailable <- IDMaleDMatureAvailable[!IDMaleDMatureAvailable %in% relatedMales] # remove them from the potential males to become breeder
      }
      
      # Select randomly one male among the selected males
      IDNewAlphaM <- sample.vec(IDMaleDMatureAvailable, 1, replace = FALSE)
      # Replace the missing breeding male by the selected one
      wolves <- NLset(turtles = wolves,
                      agents = turtle(turtles = wolves, who = IDNewAlphaM),
                      var = c("alpha", "packID", "disp"), val = cbind(alpha = 1, packID = eachPack, disp = 0))
      if(length(IDNewAlphaM) == 1){ # update the counter
        totalReplace <- totalReplace + 1
      }
    }
  }
  
  if(runTests){
    # Check if the dispersers that became breeders are no longer dispersers
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
    newNumAlpha <- NLcount(NLwith(agents = wolves, var = "alpha", val = 1))
    expect_equal(numDisp - totalReplace, newNumDisp)
    expect_equal(numAlpha + totalReplace, newNumAlpha)
  }
  
  return(wolves)
}
#######################################

########################
# Establishment in pairs
establishPairing <- function(wolves, allPackID, allWolvesRelatedness){

  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers

  # Keep in memory for the tests later
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
  wolvesPackIDStart <- wolvesPackID 
  totalEstabypairing <- 0
  
  # If the environment is not at carrying capacity, some new packs can be created
  if(CarryingCapacity > length(wolvesPackID)){

    Disperser <- NLwith(agents = wolves, var = "disp", val = 1)
    MatureDisp <- NLwith(agents = Disperser, var = "age", val = c(2:15))

    # Dispersers according to sex
    FemalesWhoEst <- NLwith(agents = MatureDisp, var = "sex", val = "F") # mature female dispersers
    IDFemaleWhoEst <- of(agents = FemalesWhoEst, var = "who")
    MalesWhoEst <- NLwith(agents = MatureDisp, var = "sex", val = "M") # matue male dispersers
    IDMaleWhoEst <- of(agents = MalesWhoEst, var = "who")
    
    # Number of possible new packs
    NumPossPacks <- CarryingCapacity - length(wolvesPackID)
    
    # Form the pairs from the female point of view (female or male does not matter as two dispersers are needed)
    # Select female dispersers who will establish using the density dependent probability of establishment
    WhoFemaleEstablished <- rbinom(n = length(IDFemaleWhoEst), size = 1, prob = NumPossPacks / CarryingCapacity) 
    IDMatDispFemWhoEstablished <- IDFemaleWhoEst[WhoFemaleEstablished == 1] # select the ID of the female dispersers who will establish
    # If there are too many females, remove the ones who cannot establish because we reach the carrying capacity
    if(length(IDMatDispFemWhoEstablished) > NumPossPacks){ 
      IDMatDispFemWhoEstablished <- sample.vec(IDMatDispFemWhoEstablished, NumPossPacks, replace = FALSE)
    }
    # Shuffle IDMatDispFemWhoEstablished so that it's not always the females with the smallest IDs that will pair first
    IDMatDispFemWhoEstablished <- sample.vec(IDMatDispFemWhoEstablished)
    
    # For each dispersing female who can establish in pair, find a dispersing male
    for(eachIndividual in IDMatDispFemWhoEstablished){
      
      # Remove from the IDMaleWhoEst the males that are too closely related
      if(length(IDMaleWhoEst) != 0){
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == eachIndividual,
                                                 colnames(allWolvesRelatedness) %in% IDMaleWhoEst, drop = FALSE]
        relatedMales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of males too closely related
        IDMaleWhoEst <- IDMaleWhoEst[!IDMaleWhoEst %in% relatedMales] # remove them from the potential males to become breeder
      }
      # Select randomly one male among the selected males
      male_disp_partner <- sample.vec(IDMaleWhoEst, 1, replace = FALSE) 
      
      # If there is a dispersing male available to pair with the dispersing female
      if(length(male_disp_partner) != 0){
        IDMaleWhoEst <- IDMaleWhoEst[!IDMaleWhoEst %in% male_disp_partner] # for next loops, remove the selected male from the possible male partners
        
        # Update the status of the two dispersers establishing together in pair
        wolves <- NLset(turtles = wolves,
                        agents = turtle(turtles = wolves, who = c(male_disp_partner, eachIndividual)), 
                        var = c("alpha", "packID", "disp"), val = cbind(alpha = rep(1, 2),
                                                                        packID = rep(max(allPackID) + 1, 2),
                                                                        disp = rep(0, 2)))
        
        # Update the counters
        allPackID <- c(allPackID, max(allPackID) + 1)
        totalEstabypairing <- totalEstabypairing + 2
      }
    }
  }
  
  if(runTests){
    # There should not be more packs than the pack carrying capacity
    newWolvesPackID <- unique(of(agents = wolves, var = "packID"))
    newWolvesPackID <- newWolvesPackID[!is.na(newWolvesPackID)] # remove the NA from the dispersers
    expect_true(length(newWolvesPackID) <= CarryingCapacity)
    # packID are always increasing unless there were no pack anymore
    if(length(wolvesPackIDStart) != 0){
      expect_equal(min(newWolvesPackID), min(wolvesPackIDStart))
    }
    expect_equal(length(allPackID), length(unique(allPackID)))
    # Dispersers that established in pairs should not be dispersers anymore
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
    expect_equal(numDisp - totalEstabypairing, newNumDisp)
    expect_equal(length(wolvesPackIDStart) + totalEstabypairing / 2, length(newWolvesPackID))
  }
  
  return(list(wolves, allPackID))
}
########################

##########################
# Establishment by budding
establishBudding <- function(wolves, allPackID, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  # Keep in memory for the tests later
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
  wolvesPackIDStart <- wolvesPackID
  totalEsta <- 0
  
  # If the environment is not at carrying capacity, some new packs can be created
  if(CarryingCapacity > length(wolvesPackID)){
    
    Disperser <- NLwith(agents = wolves, var = "disp", val = 1)
    MatureDisp <- NLwith(agents = Disperser, var = "age", val = c(2:15))
    IDMatureDisp <- of(agents = MatureDisp, var = "who")
    
    # Number of possible new packs
    NumPossPacks <- CarryingCapacity - length(wolvesPackID)
    
    # Select dispersers who will establish using the density dependent probability of establishment reduced by the probability of finding a partner to budd with (probBudd)
    WhoEstablished <- rbinom(n = length(IDMatureDisp), size = 1, prob = NumPossPacks / CarryingCapacity) & rbinom(n = length(IDMatureDisp), size = 1, prob = probBudd)
    IDMatDispWhoEstablished <- IDMatureDisp[WhoEstablished == 1] # select the ID of dispersers who will establish
    # If there are too many dispersers, remove the ones who cannot establish because we reach the carrying capacity
    if(length(IDMatDispWhoEstablished) > NumPossPacks){
      IDMatDispWhoEstablished <- sample.vec(IDMatDispWhoEstablished, NumPossPacks, replace = FALSE)
    }
    
    # Dispersers according to sex
    FemalesWhoEst <- NLwith(agents = wolves[wolves$sex == "F"], var = "who", val = IDMatDispWhoEstablished) # mature female dispersers
    IDFemaleWhoEst <- of(agents = FemalesWhoEst, var = "who")
    # Shuffle IDFemaleWhoEst so that it's not always the females with the smallest IDs that will budd first
    IDFemaleWhoEst <- sample.vec(IDFemaleWhoEst)
    MalesWhoEst <- NLwith(agents = wolves[wolves$sex == "M"], var = "who", val = IDMatDispWhoEstablished) # mature male dispersers
    IDMaleWhoEst <- of(agents = MalesWhoEst, var = "who")
    # Shuffle IDMaleWhoEst so that it's not always the males with the smallest IDs that will budd first
    IDMaleWhoEst <- sample.vec(IDMaleWhoEst)
    
    # Mature subordinate partners of the opposite sex
    FemSubPotentialMatch <- NLwith(agents = wolves[wolves$disp == 0 & wolves$alpha == 0 & wolves$age %in% 2:15,], var = "sex", val = "F") # mature female subordinates
    IDFemSubPotentialMatch <- of(agents = FemSubPotentialMatch, var = "who")
    MaleSubPotentialMatch <- NLwith(agents = wolves[wolves$disp == 0 & wolves$alpha == 0 & wolves$age %in% 2:15,], var = "sex", val = "M") # mature male subordinates
    IDMaleSubPotentialMatch <- of(agents = MaleSubPotentialMatch, var = "who")
    
    # Female disperser budding
    for(eachIndividual in IDFemaleWhoEst){
      
      # Remove from IDMaleSubPotentialMatch the males that are too closely related
      if(length(IDMaleSubPotentialMatch) != 0){
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == eachIndividual,
                                                 colnames(allWolvesRelatedness) %in% IDMaleSubPotentialMatch, drop = FALSE]
        relatedMales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of males too closely related
        IDMaleSubPotentialMatch <- IDMaleSubPotentialMatch[!IDMaleSubPotentialMatch %in% relatedMales] # remove them from the potential males to become breeder
      }
      # Select randomly one male among the selected males
      male_sub_partner <- sample.vec(IDMaleSubPotentialMatch, 1, replace = FALSE) 
      
      # If there is an available mature subordinate male to bud with the dispersing female
      if(length(male_sub_partner) != 0){
        IDMaleSubPotentialMatch <- IDMaleSubPotentialMatch[!IDMaleSubPotentialMatch %in% male_sub_partner] # for next loops, remove the selected subordinate from the possible subordinates partners
        
        # Update the status of the dispersersing female and the subordinate male who establish together by budding
        wolves <- NLset(turtles = wolves,
                        agents = turtle(turtles = wolves, who = c(male_sub_partner, eachIndividual)),
                        var = c("alpha", "packID", "disp"), val = cbind(alpha = rep(1, 2),
                                                                        packID = rep(max(allPackID) + 1, 2),
                                                                        disp = rep(0, 2)))
        
        # Update the counters
        allPackID <- c(allPackID, max(allPackID) + 1)
        totalEsta <- totalEsta + 1
      }
    }
    
    # Male disperser budding
    for(eachIndividual in IDMaleWhoEst){
      
      # Remove from IDFemSubPotentialMatch the females that are too closely related
      if(length(IDFemSubPotentialMatch) != 0){
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == eachIndividual,
                                                 colnames(allWolvesRelatedness) %in% IDFemSubPotentialMatch, drop = FALSE]
        relatedFemales <- as.numeric(colnames(relatAlphaCouple)[relatAlphaCouple[1, ] > thresholdRelatedness]) # "who" number of females too closely related
        IDFemSubPotentialMatch <- IDFemSubPotentialMatch[!IDFemSubPotentialMatch %in% relatedFemales] # remove them from the potential females to become breeder
      }
      # Select randomly one female among the selected females
      female_sub_partner <- sample.vec(IDFemSubPotentialMatch, 1, replace = FALSE) 
      
      # If there is an available mature subordinate female to bud with the dispersing male
      if(length(female_sub_partner) != 0){
        IDFemSubPotentialMatch <- IDFemSubPotentialMatch[!IDFemSubPotentialMatch %in% female_sub_partner] # for next loops, remove the selected subordinate from the possible subordinates partners
        
        # Update the status of the dispersersing male and the subordinate female who establish together by budding
        wolves <- NLset(turtles=wolves,
                        agents = turtle(turtles = wolves, who = c(female_sub_partner, eachIndividual)),
                        var = c("alpha", "packID", "disp"), val = cbind(alpha = rep(1, 2),
                                                                        packID = rep(max(allPackID) + 1, 2),
                                                                        disp = rep(0, 2)))
        
        # Update the counters
        allPackID <- c(allPackID, max(allPackID) + 1)
        totalEsta <- totalEsta + 1
      }
    }
  }
  
  if(runTests){
    # There should not be more packs than the pack carrying capacity
    newWolvesPackID <- unique(of(agents = wolves, var = "packID"))
    newWolvesPackID <- newWolvesPackID[!is.na(newWolvesPackID)] # remove the NA from the dispersers
    expect_true(length(newWolvesPackID) <= CarryingCapacity)
    # packID are always increasing unless there were no pack anymore
    if(length(wolvesPackIDStart) != 0){
      expect_equal(min(newWolvesPackID),min(wolvesPackIDStart))
    }
    expect_equal(length(allPackID), length(unique(allPackID)))
    # Dispersers that established by budding should not be dispersers anymore
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
    expect_equal(numDisp - totalEsta, newNumDisp)
    expect_equal(length(wolvesPackIDStart) + totalEsta, length(newWolvesPackID))
  }
  
  return(list(wolves, allPackID))
}
##########################

#####################
# Establishment alone
establishAlone <- function(wolves, allPackID){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  # Keep in memory for the tests later
  numDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
  wolvesPackIDStart <- wolvesPackID
  totalEsta <- 0
  
  # If the environment is not at carrying capacity, some new packs can be created
  if(CarryingCapacity > length(wolvesPackID)){
    
    # Dispersers than can establish alone
    Disperser <- NLwith(agents = wolves, var = "disp", val = 1)
    MatureDisp <- NLwith(agents = Disperser, var = "age", val = c(2:15))
    IDMatureDisp <- of(agents = MatureDisp, var = "who")
    
    # Number of possible new packs that can be created
    NumPossPacks <- CarryingCapacity - length(wolvesPackID)
    
    # Select dispersers who will establish using the density dependent probability of establishment
    WhoEstablished <- rbinom(n = length(IDMatureDisp), size = 1, prob = NumPossPacks / CarryingCapacity) 
    IDMatDispWhoEstablished <- IDMatureDisp[WhoEstablished == 1] # select the ID of dispersers who will establish
    # If there are too many dispersers, remove the ones who cannot establish because we reach the carrying capacity
    if(length(IDMatDispWhoEstablished) > NumPossPacks){
      IDMatDispWhoEstablished <- sample.vec(IDMatDispWhoEstablished, NumPossPacks, replace = FALSE)
    }
    
    for(eachIndividual in IDMatDispWhoEstablished){ 
      
      # Update the status of the dispersers that establish alone
      wolves <- NLset(turtles = wolves, 
                      agents = turtle(turtles = wolves, who = eachIndividual),
                      var = c("alpha", "packID", "disp"), val = cbind(alpha = 1,
                                                                      packID = max(allPackID) + 1,
                                                                      disp = 0))
      
      # Update counters
      allPackID <- c(allPackID, max(allPackID) + 1)
      totalEsta <- totalEsta + 1
    }
  }
  
  if(runTests){
    # There should not be more packs than the pack carrying capacity
    newWolvesPackID <- unique(of(agents = wolves, var = "packID"))
    newWolvesPackID <- newWolvesPackID[!is.na(newWolvesPackID)] # remove the NA from the dispersers
    expect_true(length(newWolvesPackID) <= CarryingCapacity)
    # packID are always increasing unless there were no pack anymore
    if(length(wolvesPackIDStart) != 0){
      expect_equal(min(newWolvesPackID),min(wolvesPackIDStart))
    }
    expect_equal(length(allPackID), length(unique(allPackID)))
    # Dispersers that established alone should not be dispersers anymore
    newNumDisp <- NLcount(NLwith(agents = wolves, var = "disp", val = 1))
    expect_equal(numDisp - totalEsta, newNumDisp)
    expect_equal(length(wolvesPackIDStart) + totalEsta, length(newWolvesPackID))
  }
  
  return(list(wolves, allPackID))
}
#####################

###############################################
# Replacement of breeding males by subordinates
MaleAlphaSurbordinateReplacement <- function(wolves, allWolvesRelatedness){
  
  wolvesPackID <- unique(of(agents = wolves, var = "packID"))
  wolvesPackID <- wolvesPackID[!is.na(wolvesPackID)] # remove the NA from the dispersers
  
  for(eachPack in wolvesPackID){
    
    wolvesInPack <- NLwith(agents = wolves, var = "packID", val = eachPack)
    MaleInPack <- NLwith(agents = wolvesInPack, var = "sex", val = "M") # males in the pack
    MaleAlphaInPack <- NLwith(agents = MaleInPack, var = "alpha", val = 1) # breeding male in the pack
    MaleSMatureAvailable <- NLwith(agents = MaleInPack, var = "age", val = c(2:15)) # mature subordinates in the pack
    
    # If there is no breeding male but other males available in the pack
    if(NLcount(MaleAlphaInPack) == 0 & NLcount(MaleSMatureAvailable) != 0){
      IDMaleSMatureAvailable <- of(agents = MaleSMatureAvailable, var = "who")
      
      # If there are several males to replace the missing breeding male, take the one the least related to the current breeding female
      femaleAlphaInPack <- NLwith(agents = wolvesInPack, var = "alpha", val = 1) # identify the breeding female
      if(NLcount(MaleSMatureAvailable) > 1 & NLcount(femaleAlphaInPack) == 1){
        IDfemaleAlphaInPack <- of(agents = femaleAlphaInPack, var = "who")
        # Relatedness between the breeders
        relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDfemaleAlphaInPack,
                                                 colnames(allWolvesRelatedness) %in% IDMaleSMatureAvailable]
        lessRelated <- which(relatAlphaCouple == min(relatAlphaCouple)) # male(s) least related to the current alpha female
        
        # If there are several males with the same relatedness
        if(length(lessRelated) > 1){
          
          # First, check if there is among them a dismissed male, and if yes, take this one back as breeding
          if(sum(of(agents = turtle(turtles = wolves, who = IDMaleSMatureAvailable[lessRelated]), var = "dismissed")) == 1){
            IDNewAlphaM <- of(agents = NLwith(agents = wolvesInPack, var = "dismissed", val = 1), var = "who")
          } else {
            # If there was no dismissed breeding male, choose one among the least related at random
            lessRelated <- sample.vec(lessRelated, 1)
            IDNewAlphaM <- IDMaleSMatureAvailable[lessRelated]
          }
        } else {
          # If there is only one male least related, take this one
          IDNewAlphaM <- IDMaleSMatureAvailable[lessRelated]
        }
      } else {
        # If there is no breeding female, select randomly one male among the selected males
        IDNewAlphaM <- sample.vec(IDMaleSMatureAvailable, 1, replace = FALSE) 
        
        if(runTests){
          # If there is no female in the pack, there should be no any dismissed breeding male in the pack
          if(NLcount(femaleAlphaInPack) == 0){
            expect_equivalent(NLcount(NLwith(agents = wolvesInPack, var = "dismissed", val = 1)), 0)
          }
        }
        
      }
      
      # Update the status of the new breeding male
      wolves <- NLset(turtles = wolves,
                      agents = turtle(turtles = wolves, who = IDNewAlphaM),
                      var = "alpha", val = 1)
      
      # Check the relatedness between the two breeding individuals (if the male was not a previously dismissed)
      if(of(agents = turtle(turtles = wolves, who = IDNewAlphaM), var = "dismissed") != 1){
        # Calculate the relatedness between the two alpha
        if(NLcount(femaleAlphaInPack) == 1){
          IDfemaleAlphaInPack <- of(agents = femaleAlphaInPack, var = "who")
          # Relatedness between the breeders
          relatAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) == IDfemaleAlphaInPack,
                                                   colnames(allWolvesRelatedness) == IDNewAlphaM]
          
          # If the breeders are too related, check if there is not a female in the pack less related to the new breeding male
          if(relatAlphaCouple > thresholdRelatedness){ 
            allFemalesInPack <- NLwith(agents = wolvesInPack, var = "sex", val = "F") # females in the pack
            IDallFemalesInPack <- of(agents = NLwith(agents = allFemalesInPack, var = "age", val = 2:16), var = "who") # mature females in the pack
            
            # If there are mature females other than the current breeding female
            if(length(IDallFemalesInPack) > 1){ 
              relatPotentialAlphaCouple <- allWolvesRelatedness[rownames(allWolvesRelatedness) %in% IDallFemalesInPack,
                                                                colnames(allWolvesRelatedness) == IDNewAlphaM]
              
              # If the current couple is more related than other potential couples
              if(sum(relatAlphaCouple > relatPotentialAlphaCouple) != 0){
                IDnewAlphaF <- as.numeric(names(which(relatPotentialAlphaCouple == min(relatPotentialAlphaCouple)))) # id of the least related female
                if(length(IDnewAlphaF) > 1){ # if there are several least related females
                  IDnewAlphaF <- sample.vec(IDnewAlphaF, 1) # take one randomly
                }
                
                # Make the least related female, the new breeding female and put the former breeding female as a subordinate
                wolves <- NLset(turtles = wolves, agents = turtle(turtles = wolves, who = IDnewAlphaF),
                                var = "alpha", val = 1)
                wolves <- NLset(turtles = wolves, agents = femaleAlphaInPack,
                                var = "alpha", val = 0)
                
                if(runTests){
                  # The new female should not be the same as the former one
                  expect_false(IDnewAlphaF == of(agents = femaleAlphaInPack, var = "who"))
                }
                
              }
            }
          }
        }
      }
    }
  }
  
  # Remove all infos about dismissed breeding position as all replacement have been done
  wolves <- NLset(turtles = wolves, agents = wolves, var = "dismissed", val = 0)
  
  if(runTests){
    # Any packs with mature males (>= 2 yrs old) must have a breeding male
    # There must be as many unique packID with mature males in it, than the number of breeding males in these packs
    allMale <- NLwith(agents = wolves, var = "sex", val = "M")
    allMaleInPack <- NLwith(agents = allMale, var = "disp", val = 0)
    matureMaleInPack <- NLwith(agents = allMaleInPack, var = "age", val = 2:15)
    expect_equal(length(unique(of(agents = matureMaleInPack, var = "packID"))),
                 NLcount(NLwith(agents = matureMaleInPack, var = "alpha", val = 1)))
    # There should be no dismissed individuals anymore
    expect_equivalent(NLcount(NLwith(agents = wolves, var = "dismissed", val = 1)), 0)
    # There should be only one breeding female and one breeding male per pack
    males <- NLwith(agents = wolves, var = "sex", val = "M")
    females <- NLwith(agents = wolves, var = "sex", val = "F")
    maleAlpha <- NLwith(agents = males, var = "alpha", val = 1)
    femaleAlpha <- NLwith(agents = females, var = "alpha", val = 1)
    expect_equivalent(length(unique(of(agents = maleAlpha, var = "packID"))), NLcount(maleAlpha))
    expect_equivalent(length(unique(of(agents = femaleAlpha, var = "packID"))), NLcount(femaleAlpha))
  }
  
  return(wolves)
}
###############################################
