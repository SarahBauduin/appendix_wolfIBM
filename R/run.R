########################
# Wolf IBM, May 2020
# Model authors (in alphabetical order): Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi
#
# From the submitted paper: Exploring the impact of lesser-known social dynamics on wolf populations through an individual-based approach
# Paper authors: Sarah Bauduin, Oksana Grente, Nina Luisa Santostasi, Paolo Ciucci, Christophe Duchamp, and Olivier Gimenez
# Submitted to Ecological Modelling
########################


########################
# Run the complete model
########################


###################
# Prepare the model

# Call the sub-models
source("submodels.R") # the functions in this file do not need to be modified
# The package pedantics is supposed to be installed and loaded now
# If it did not work, please verify the path to the package folder (/appendix_wolfIBM/R/pedantics_1.7.tar.gz) in the file submodels.R

# Call the function to init the population and the sub-model parameter values
source("initParam.R") # the data about the initial population and the sub-model parameters may be changed to be adapted to the user's question

# Define the model parameter
nReplicate <- 200 # how many replicates of the population simulation
nYearSim <- 25 # how many years of simulation

# Create the output files
popSimRep <- list() # record as the list the state of the population each year for each simulation
# Create a df to record the dynamic of the packs (i.e., creation and losses)
packDyn <- data.frame(repSim = rep(1:nReplicate, each = nYearSim),
                      yearSim = rep(1:nYearSim, nReplicate),
                      numPack = 0, # number of packs at the beginning of the year
                      lostPackMort = 0, # number of packs that disappeared during the year because all the individuals in the pack died 
                      packDiss = 0, # number of packs that dissolved during the year
                      newPackPair = 0, # number of packs created during the year by the pairing of two dispersing individuals
                      newPackBud = 0, # number of packs created during the year by "budding"
                      newPackSingle = 0) # number of packs created during the the year by individuals alone
###################


###########################
# Run the model with a 
# Model version M1

for(j in 1:nReplicate){ # loop over the simulation replicates
  
  # Create the initial population at the beginning of each simulation replicate
  wolves <- init(initPopWolf_DF)
  # Create a list to record the state of the population each year (for THIS replicate)
  popSim <- list()
  popSim[[1]] <- wolves # record the initial state of the wolf population
  
  # Create some vectors needed for the sub-models
  allWolvesID <- of(agents = wolves, var = "who") # keep in memory all the wolves ID ever created
  allPackID <- unique(of(agents = wolves, var = "packID")) # keep in memory all the pack ID ever created
  allPackID <- allPackID[!is.na(allPackID)] # remove NA (i.e., packID for dispersers)
  yearSim <- 0 # used to indicate the cohorts, the number itself doesn't matter, it just needs to be updated at each loop
  
  for(i in 1:nYearSim){ # loop over the years simulated
    
    # Update the year simulated at the beginning of each year
    yearSim <- yearSim + 1
    
    # Record the pack dynamics
    packIDnow1 <- unique(of(agents = wolves, var = "packID"))
    packIDnow1 <- packIDnow1[!is.na(packIDnow1)]
    packDyn[packDyn$repSim == j & packDyn$yearSim == i, "numPack"] <- length(packIDnow1) # how many packs there is at the beginning of this year
    
    if(NLcount(wolves) != 0){ # do not execute the sub-models if there are no wolves alive
      
      ##############
      # Reproduction
      resRepro <- repro(wolves, allWolvesID, yearSim, popSim) # result = list(wolves, allWolvesID, allWolvesRelatedness). yearSim  and popSim are used but not modified so not returned
      wolves <- resRepro[[1]] # update the objects each time there are modified
      allWolvesID <- resRepro[[2]]
      allWolvesRelatedness <- resRepro[[3]]
      ##############
      
      #######
      # Aging
      wolves <- aging(wolves)
      #######
      
      ###########
      # Mortality
      wolves <- mortality(wolves)
      # Record pack dynamics after mortality as some packs may have disappeared if all wolves from a pack died
      packIDnow2 <- unique(of(agents = wolves, var = "packID"))
      packIDnow2 <- packIDnow2[!is.na(packIDnow2)]
      packDyn[packDyn$repSim == j & packDyn$yearSim == i, "lostPackMort"] <- length(packIDnow1) - length(packIDnow2) # how many packs were lost after the mortality process
      ###########
    }
    
    if(NLcount(wolves) != 0){ # do not execute the sub-models if there are no wolves alive (i.e., if they all died in the mortality sub-model)
      
      ###################
      # Pack dissolvement
      wolves <- PackDissolvement(wolves)
      # Record pack dynamics after pack dissolvement as some packs may have disappeared
      packIDnow3 <- unique(of(agents = wolves, var = "packID"))
      packIDnow3 <- packIDnow3[!is.na(packIDnow3)]
      packDyn[packDyn$repSim == j & packDyn$yearSim == i, "packDiss"] <- length(packIDnow2) - length(packIDnow3) # how many packs were lost after pack dissolvement
      ###################
      
      #################################################
      # Replacement of breeding females by subordinates
      wolves <- FemaleAlphaSurbordinateReplacement(wolves, allWolvesRelatedness)
      #################################################
      
      ###########
      # Dispersal
      resDispersal <- dispersal(wolves) # result = list(wolves, max_sizes)
      wolves <- resDispersal[[1]]
      max_sizes <- resDispersal[[2]]
      ###########
      
      ########################
      # Immigration/Emigration
      resImmigr <- immigration(wolves, allWolvesID, popSim) # result = list(wolves, allWolvesID, allWolvesRelatedness)
      wolves <- resImmigr[[1]]
      allWolvesID <- resImmigr[[2]]
      allWolvesRelatedness <- resImmigr[[3]]
      wolves <- emigration(wolves)
      ########################
      
      ##########
      # Adoption
      wolves <- Adoptee(wolves, max_sizes)
      ##########
      
      #######################################
      # Replacement of breeders by dispersers
      wolves <- AlphaDisperserReplacement(wolves, allWolvesRelatedness)
      #######################################
      
      ########################
      # Establishment in pairs
      resEstaPairing <- establishPairing(wolves, allPackID, allWolvesRelatedness) # result = list(wolves, allPackID)
      wolves <- resEstaPairing[[1]]
      allPackID <- resEstaPairing[[2]]
      # Record pack dynamics after the establishment in pairs as some packs may have been created
      packIDnow4 <- unique(of(agents = wolves, var = "packID"))
      packIDnow4 <- packIDnow4[!is.na(packIDnow4)]
      packDyn[packDyn$repSim == j & packDyn$yearSim == i, "newPackPair"] <- length(packIDnow4) - length(packIDnow3) # how many packs were added after the establishment in pairs
      ########################
      
      ##########################
      # Establishment by budding
      resEstaBudding <- establishBudding(wolves, allPackID, allWolvesRelatedness) # result = list(wolves, allPackID)
      wolves <- resEstaBudding[[1]]
      allPackID <- resEstaBudding[[2]] 
      # Record pack dynamics after the establishment by budding as some packs may have been created
      packIDnow5 <- unique(of(agents = wolves, var = "packID"))
      packIDnow5 <- packIDnow5[!is.na(packIDnow5)]
      packDyn[packDyn$repSim == j & packDyn$yearSim == i, "newPackBud"] <- length(packIDnow5) - length(packIDnow4) # how many packs were added after the establishment by budding
      ##########################
      
      #####################
      # Establishment alone
      resEstaAlone <- establishAlone(wolves, allPackID) # result = list(wolves, allPackID)
      wolves <- resEstaAlone[[1]]
      allPackID <- resEstaAlone[[2]]
      # Record pack dynamics  after the establishment alone as some packs may have been created
      packIDnow6 <- unique(of(agents = wolves, var = "packID"))
      packIDnow6 <- packIDnow6[!is.na(packIDnow6)]
      packDyn[packDyn$repSim == j & packDyn$yearSim == i, "newPackSingle"] <- length(packIDnow6) - length(packIDnow5) # how many packs were added after the establishment alone
      #####################
      
      ###############################################
      # Replacement of breeding males by subordinates
      wolves <- MaleAlphaSurbordinateReplacement(wolves, allWolvesRelatedness)
      ###############################################
      
      # Add the new state of the population
      popSim[[i + 1]] <- wolves
      
    } else { # if there are no wolves anymore
      
      # Do not run the sub-models but still add the new state of the population (i.e., no wolves)
      popSim[[i + 1]] <- wolves
    }
    
  } # end of the year simulated
  
  # Add the list popSim to the list popSimRep
  popSimRep[[j]] <- popSim # population simulated over all the years for the current simulation replicate
  # Print the number of the current replicate
  print(paste0("Replicate", j))
  
  # At the end of each replicate, save popSimRep and packDyn in a .Rdata file so that if there a problem, you can access to all simulation replicates already done
  save(popSimRep, packDyn, file = "wolfSimRes.RData")
  
} # end of the replicate simulation replicate
###########################


############################
# Exploration of the results
load("wolfSimRes.RData")

#######################################
# Result trend over the simulated years
resSim_pack2breeders <- list() # number of packs with 2 breeding individuals, length(resSim_pack2breeders) == nYearSim
for(i in 1:(nYearSim + 1)){ # nYearSim + 1 because the population was recorded before the beginning of the simulation and the population was recorded at the end of each year simulated
  j <- 1
  # Number of packs that have 2 breeding individuals in them for this year simulated for this simulation replicate
  resSim_pack2breeders[[i]] <- length(which(table(of(agents = NLwith(agents = popSimRep[[j]][[i]], var = "alpha", val = 1), var = "packID")) == 2)) 
  for(j in 2:nReplicate){
    resSim_pack2breeders[[i]] <- c(resSim_pack2breeders[[i]], length(which(table(of(agents = NLwith(agents = popSimRep[[j]][[i]], var = "alpha", val = 1), var = "packID")) == 2)))
  }
}

# Plot the mean number of packs with a breeding pair and the 95% confidence interval around the mean
plot(0:nYearSim, unlist(lapply(resSim_pack2breeders, mean)), type = "l",
     xlab = "Years simulated", ylab = "Number of packs with a breeding pair",
     ylim = c(min(unlist(resSim_pack2breeders)), max(unlist(resSim_pack2breeders))))
lines(0:nYearSim, unlist(lapply(resSim_pack2breeders, quantile, probs = 0.025)), lty = 2) # dashed line
lines(0:nYearSim, unlist(lapply(resSim_pack2breeders, quantile, probs = 0.975)), lty = 2) # dashed line

#######################################

#########################################
# Results for the last year of simulation
resSim_nInd <- numeric() # number of individuals
resSim_newPacks <- numeric() # number of new packs created
resSim_propRes <- numeric() # proportion of resident individuals
resSim_relatBreeders <- numeric() # relatedness between breeders in packs

for(j in 1:nReplicate){
  # Total number of individuals for the last year of simulation for this simulation replicate
  resSim_nInd <- c(resSim_nInd, NLcount(popSimRep[[j]][[(nYearSim + 1)]]))
  # Sum of new packs created by dispersing individuals pairing, by budding and by dispersers alone for the last year of simulation for this simulation replicate
  resSim_newPacks <- c(resSim_newPacks, sum(packDyn[packDyn$repSim == j & packDyn$yearSim == nYearSim, c("newPackPair", "newPackBud",  "newPackBud")]))
  # Number of resident individuals (i.e., not dispersers) over the total number of individual for the last year of simulation for this simulation replicate
  resSim_propRes <- c(resSim_propRes, NLcount(NLwith(agents = popSimRep[[j]][[(nYearSim + 1)]], var = "disp", val = 0)) / NLcount(popSimRep[[j]][[(nYearSim + 1)]]))
  
  # Relatedness between breeders in packs for the last year of simulation for this simulation replicate
  # Use of the function relatedness() available in submodels.R, created using functions from the package pedantics
  relatednessCouples <- NA
  withParentID <- NLwith(agents = popSimRep[[j]][[(nYearSim + 1)]], var = "motherID", val = 0:1000000) # extract individuals with parent ID
  alpha <- NLwith(agents = withParentID, var = "alpha", val = 1)
  alphaFemale <- NLwith(agents = alpha, var = "sex", val = "F") # female breeders
  alphaMale <- NLwith(agents = alpha, var = "sex", val = "M") # male breeders
  packID2alphas <- intersect(of(agents = alphaFemale, var = "packID"), of(agents = alphaMale, var = "packID")) # packs with a breeding pair
  allRelatedness <- relatedness(listAllInd = popSimRep[[j]], whoInd = of(agents = alpha, var = "who"))
  for(eachPack in packID2alphas){
    # Extract for each pack with a breeding pair, the relatedness between the breeding male and the breeding female of the pack
    relatedCouplePack <- allRelatedness[rownames(allRelatedness) == of(agents = NLwith(agents = alphaFemale, var = "packID", val = eachPack), var = "who"),
                                        colnames(allRelatedness) == of(agents = NLwith(agents = alphaMale, var = "packID", val = eachPack), var = "who")]
    relatednessCouples <- c(relatednessCouples, relatedCouplePack)
  }
  resSim_relatBreeders <- c(resSim_relatBreeders, relatednessCouples)
}

# Plot a summary of the results for the last year of simation over all simulation replicates using boxplot
boxplot(resSim_nInd, main = "Number of individuals in the population")
boxplot(resSim_newPacks, main = "Number of new packs formed in the year")
boxplot(resSim_propRes, main = "Proportion of resident individuals in the population")
boxplot(resSim_relatBreeders, main = "Relatedness between the breeders in the packs")
#########################################

############################
