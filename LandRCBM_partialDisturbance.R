defineModule(sim, list(
  name = "LandRCBM_partialDisturbance",
  description = "",
  keywords = "",
  authors = structure(list(list(given = "Camille", family = "Giuliano", role = c("aut", "cre"), email = "cams0405@live.ca", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(LandRCBM_partialDisturbance = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "LandRCBM_partialDisturbance.Rmd"),
  reqdPkgs = list("SpaDES.core (>= 2.1.8.9018)", "data.table", "PredictiveEcology/LandR@development"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                    "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput("cohortData", "data.table", 
                 desc = "age cohort-biomass table hooked to pixel group map by pixelGroupIndex at succession time step"),
    expectsInput("PartialDistTable", "data.table",
                 desc ="Table with partial disturbance information"),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "updated community map at each succession time step"),
    expectsInput("species", "data.table",
                 desc = "a table that has species traits such as longevity...",
                 sourceURL = "https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/species.txt"),
    expectsInput("speciesEcoregion", "data.table",
                 desc = "table defining the maxANPP, maxB and SEP, which can change with both ecoregion and simulation time",
                 sourceURL = "https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession-dynamic-inputs_test.txt"),
    # expectsInput("sufficientLight", "data.frame",
    #              desc = "table defining how the species with different shade tolerance respond to stand shadeness",
    #              sourceURL = "https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession_test.txt"),
    # expectsInput("treedFirePixelTableSinceLastDisp", "data.table",
    #              desc = "3 columns: pixelIndex, pixelGroup, and burnTime. Each row represents a forested pixel that was burned up to and including this year, since last dispersal event, with its corresponding pixelGroup and time it occurred")
  ),
  outputObjects = bindrows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput("cohortData", "data.table",
                  desc = paste("age cohort-biomass table hooked to pixel group map",
                               "by pixelGroupIndex at succession time step")),
    createsOutput("pixelGroupMap", "RasterLayer",
                  desc = "updated community map at each succession time step"),
    # createsOutput("serotinyResproutSuccessPixels", "numeric",
    #               desc = "Pixels that were successfully regenerated via serotiny or resprouting. This is a subset of treedBurnLoci")
  )
))

doEvent.LandRCBM_partialDisturbance = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim), "LandRCBM_partialDisturbance", "annualPartialDist", eventPriority = 3)
    },
    
    annualPartialDist = {
      sim <- processDist(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "LandRCBM_partialDisturbance", "annualPartialDist", eventPriority = 4)
    },
    
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  
  return(invisible(sim))
}

processDist <- function(sim) {
  if (as.numeric(time(sim)) %in% sim$partialDistTable$distYear) {
    if (is.null(sim$partialDistLoc)) {
      partialDistTable <- sim$partialDistTable
      partialDist <- processPartialDist(cohortData = sim$cohortData, 
                                        partialDistTable = as.data.table(partialDistTable), 
                                        pixelGroupMap = sim$pixelGroupMap,
                                        currentTime = time(sim))
    } else {
      partialDistTable <- sim$partialDistTable
      partialDist <- processPartialDist(cohortData = sim$cohortData, 
                                        partialDistTable = as.data.table(partialDistTable), 
                                        pixelGroupMap = sim$pixelGroupMap,
                                        currentTime = time(sim),
                                        partialDistLoc = sim$partialDistLoc) }
    
    ## update cohortData and pixelGroupMap
    sim$cohortData <- partialDist$cohortData
    sim$pixelGroupMap <- partialDist$pixelGroupMap
    
    #eventually probably have 3 functions: 
    #harvest (partial dist with replanting), 
    #fire (partial dist with severity calcs that decide what dies, resprouting +serotiny),
    #other (partial dist with resprouting)
    
  }
  
  #if run with simpleHarvest. lots of hardcoding here, would like to make it more generic.
  if (!is.null(sim$rstCurrentHarvest) && global(sim$rstCurrentHarvest, "max", na.rm = TRUE)[[1]] > 0) {
    # browser()
    for (species in names(sim$speciesHarvestMaps)){
      
      partialDistTable <- data.table(
        speciesCode = species,
        distYear    = time(sim))
      
      partialDist <- processPartialDist(cohortData = sim$cohortData, 
                                        partialDistTable = as.data.table(partialDistTable), 
                                        pixelGroupMap = sim$pixelGroupMap,
                                        partialDistLoc = sim$speciesHarvestMaps[[species]],
                                        currentTime = time(sim))
      sim$cohortData <- partialDist$cohortData
      sim$pixelGroupMap <- partialDist$pixelGroupMap
    }
    # browser()
  }
  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  # input species ecoregion dynamics table
  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    sim$speciesEcoregion <- prepInputsSpeciesEcoregion(url = extractURL("speciesEcoregion"),
                                                       dPath = dPath, cacheTags = cacheTags)
  }
  
  
  ## ALL OBJECTS BELOW ARE FROM BIOMASS_REGENERATIONPM
  ## get LANDISII main input table where species and light requirements tables come from
  if (!suppliedElsewhere("sufficientLight", sim) |
      (!suppliedElsewhere("species", sim))) {
    mainInput <- prepInputsMainInput(url = NULL, dPath, cacheTags) ## uses default URL
  }
  
  ## read species txt and convert it to data table
  if (!suppliedElsewhere("species", sim)) {
    sim$species <- prepInputsSpecies(url = extractURL("species"), dPath, cacheTags)
  }
  
  ## make light requirements table
  if (!suppliedElsewhere("sufficientLight", sim)) {
    sufficientLight <- data.frame(mainInput)
    startRow <- which(sufficientLight$col1 == "SufficientLight")
    sufficientLight <- sufficientLight[(startRow + 1):(startRow + 5), 1:7]
    sufficientLight <- data.table(sufficientLight)
    sufficientLight <- sufficientLight[, lapply(.SD, function(x) as.numeric(x))]
    
    names(sufficientLight) <- c("speciesshadetolerance",
                                "X0", "X1", "X2", "X3", "X4", "X5")
    sim$sufficientLight <- data.frame(sufficientLight)
  }
  
  return(invisible(sim))
}