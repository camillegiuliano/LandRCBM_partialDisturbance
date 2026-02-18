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
    expectsInput("inactivePixelIndex", "logical", 
                 desc = "internal use. Keeps track of which pixels are inactive"),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "updated community map at each succession time step"),
    expectsInput("species", "data.table",
                 desc = "a table that has species traits such as longevity...",
                 sourceURL = "https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/species.txt"),
    expectsInput("speciesEcoregion", "data.table",
                 desc = "table defining the maxANPP, maxB and SEP, which can change with both ecoregion and simulation time",
                 sourceURL = "https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession-dynamic-inputs_test.txt"),
    expectsInput("sufficientLight", "data.frame",
                 desc = "table defining how the species with different shade tolerance respond to stand shadeness",
                 sourceURL = "https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession_test.txt"),
    expectsInput("treedFirePixelTableSinceLastDisp", "data.table",
                 desc = "3 columns: pixelIndex, pixelGroup, and burnTime. Each row represents a forested pixel that was burned up to and including this year, since last dispersal event, with its corresponding pixelGroup and time it occurred")
  ),
  outputObjects = bindrows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput("cohortData", "data.table",
                  desc = paste("age cohort-biomass table hooked to pixel group map",
                               "by pixelGroupIndex at succession time step")),
    createsOutput("pixelGroupMap", "RasterLayer",
                  desc = "updated community map at each succession time step"),
    createsOutput("serotinyResproutSuccessPixels", "numeric",
                  desc = "Pixels that were successfully regenerated via serotiny or resprouting. This is a subset of treedBurnLoci")
  )
))

doEvent.LandRCBM_partialDisturbance = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  
  ##TODO: checks
  browser()
  #Can bring in either a rasters or a tables of partial disturbances.
  #if raster, need to be turned into a table. Looks like the majority is actually vector data.
  
  
  if (as.numeric(time(sim)) %in% sim$partialDistTable$distYear) {
    
    
    processPartialDist(cohortData, PartialDistTable, inactivePixelIndex)
    
    
  }
  
  
  #########################
  ##EVERYTHING BELOW THIS POINT: is copy pasted from Biomass_regenerationPM
  # here site shade is calculated, resprouting occurs, new cohorts are created.
  #########################
  
  ## CALCULATE SIDE SHADE -----------------------------
  siteShade <- data.table(calcSiteShade(currentTime = round(time(sim)), burnedPixelCohortData,
                                        sim$speciesEcoregion, sim$minRelativeB))
  siteShade <- siteShade[, .(pixelGroup, siteShade)]
  burnedPixelCohortData <- siteShade[burnedPixelCohortData, on = "pixelGroup", nomatch = NA]
  burnedPixelCohortData[is.na(siteShade), siteShade := 0]
  rm(siteShade)
  
  ## clean burnedPixelCohortData from unnecessary columns
  # set(burnedPixelCohortData, NULL, c("B", "mortality", "aNPPAct"), NULL)
  set(burnedPixelCohortData, NULL, c("mortality", "aNPPAct"), NULL)
  # set(burnedPixelCohortData, ,c("sumB", "siteShade"), 0) # assume the fire burns all cohorts on site
  setkey(burnedPixelCohortData, speciesCode)
  
  ## DO RESPROUTING --------------------------
  ## assess resprouting reproduction:
  ## basically same thing as serotiny
  resproutingOutputs <- doResprouting(serotinyPixel = serotinyPixel,
                                      treedFirePixelTableSinceLastDisp = treedFirePixelTableSinceLastDisp,
                                      burnedPixelCohortData = burnedPixelCohortData,
                                      postFirePixelCohortData = postFirePixelCohortData,
                                      currentTime = time(sim), species = sim$species,
                                      sufficientLight = sim$sufficientLight,
                                      calibrate = P(sim)$calibrate,
                                      postFireRegenSummary = sim$postFireRegenSummary)
  
  postFirePixelCohortData <- resproutingOutputs$postFirePixelCohortData
  sim$serotinyResproutSuccessPixels <- resproutingOutputs$serotinyResproutSuccessPixels
  if (!is.null(resproutingOutputs$postFireRegenSummary))
    sim$postFireRegenSummary <- resproutingOutputs$postFireRegenSummary
  
  rm(resproutingOutputs)
  
  ## ADD NEW COHORTS -----------------------------
  ## add new cohorts to pixels where serotiny/regeneration were activated
  if (NROW(postFirePixelCohortData)) {
    ## redo post-fire pixel groups by adding the maxPixelGroup to their ecoregioMap values
    if (!is.null(sim$serotinyResproutSuccessPixels)) {
      
      # Add new cohorts to BOTH the sim$cohortData and sim$pixelGroupMap
      ## reclassify pixel groups as burnt (0L)
      if (verbose > 0)
        message(blue("Post serotiny and resprouting"))
      
      ## add the survivors cohorts to the serotiny/reprouting ones
      cols <- c("pixelGroup", "pixelIndex", "speciesCode", "ecoregionGroup", "age", "B")
      postFirePixelCohortData <- rbind(postFirePixelCohortData, burnedPixelCohortData[B > 0, ..cols],
                                       use.names = TRUE, fill = TRUE)
      postFirePixelCohortData[is.na(type), type := "survivor"]
      
      ## set ages to 1 here, because updateCohortData will only so so if there isn't an age column
      postFirePixelCohortData[is.na(age), age := 1L]
      
      ## filter cohortData to only have unburnt pixels -- this is not sufficient!!!
      ## in PGs where cohorts die in one but not other pixels, these cohorts from other pixels are added back where they were supposed to be removed.
      # unburnedPCohortData <- addPixels2CohortData(copy(sim$cohortData), sim$pixelGroupMap)
      # unburnedPCohortData <- unburnedPCohortData[!pixelIndex %in% treedFirePixelTableSinceLastDisp$pixelIndex]
      # set(unburnedPCohortData, NULL, "pixelIndex", NULL)  ## collapse pixel groups again
      # unburnedPCohortData <- unburnedPCohortData[!duplicated(unburnedPCohortData)]
      
      ## redo PGs in all burnt pixels
      ## 1) we need to create a table of unburt pixels and burnt pixels with dead and surviving cohorts,
      ## but not new cohorts (serotiny/resprout) -- these are added by updateCohortData
      ## 2) then remove dead cohorts for updateCohortData and redo PG
      ## the PGs need to be done twice otherwise, once to account for cohorts that only died in some but not all pixels of a given
      ## pixelGroup, and the second time to ensure that pixels that became similar after the death of some cohorts can
      ## be grouped together.
      
      unburnedPCohortData <- addPixels2CohortData(copy(sim$cohortData), sim$pixelGroupMap)
      unburnedPCohortData <- unburnedPCohortData[!pixelIndex %in% treedFirePixelTableSinceLastDisp$pixelIndex]
      newPCohortData <- rbind(unburnedPCohortData, burnedPixelCohortData, fill = TRUE)
      
      cd <- newPCohortData[, c("pixelIndex", columnsForPixelGroups), with = FALSE]
      newPCohortData[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = 0L, columns = columnsForPixelGroups)]
      pixelGroupMap <- sim$pixelGroupMap
      pixelGroupMap[newPCohortData$pixelIndex] <- newPCohortData$pixelGroup
      
      if (isTRUE(getOption("LandR.assertions", TRUE))) {
        test <- setdiff(which(!is.na(pixelGroupMap[])), newPCohortData$pixelIndex)
        if (any(pixelGroupMap[test] != 0)) {
          stop("Bug in Biomass_regenerationPM: not all pixels are in the joint burnt and unburnt pixelCohortData table")
        }
      }
      
      ## remove dead cohorts and re-do pixelGroups
      newPCohortData <- newPCohortData[B > 0]
      cd <- newPCohortData[, c("pixelIndex", columnsForPixelGroups), with = FALSE]
      newPCohortData[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = 0L, columns = columnsForPixelGroups)]
      pixelGroupMap[newPCohortData$pixelIndex] <- newPCohortData$pixelGroup
      
      ## recalculate sumB
      newPCohortData[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]
      
      ## check for duplicates at pixel-level
      if (isTRUE(getOption("LandR.assertions", TRUE))) {
        if (any(duplicated(newPCohortData[, .(speciesCode, age, pixelIndex)]))) {
          stop("Duplicate cohorts in pixels were found after burning, serotiny and resprouting")
        }
      }
      
      ## collapse to PGs
      tempCohortData <- copy(newPCohortData)
      set(tempCohortData, NULL, "pixelIndex", NULL)
      cols <- c("pixelGroup", "speciesCode", "ecoregionGroup", "age")
      tempCohortData <- tempCohortData[!duplicated(tempCohortData[, ..cols])]
      
      if (isTRUE(getOption("LandR.assertions", TRUE))) {
        if (any(duplicated(tempCohortData[, .(speciesCode, age, pixelGroup)]))) {
          stop("Duplicate cohorts in pixelGroups were found after burning, serotiny and resprouting")
        }
      }
      
      outs <- updateCohortData(newPixelCohortData = copy(postFirePixelCohortData[, -"pixelGroup", with = FALSE]),
                               cohortData = copy(tempCohortData),
                               pixelGroupMap = pixelGroupMap,
                               currentTime = round(time(sim)),
                               speciesEcoregion = copy(sim$speciesEcoregion),
                               treedFirePixelTableSinceLastDisp = copy(treedFirePixelTableSinceLastDisp),
                               initialB = P(sim)$initialB,
                               successionTimestep = P(sim)$successionTimestep)
      
      assertPostFireDist(cohortDataOrig = tempCohortData, pixelGroupMapOrig = pixelGroupMap,
                         cohortDataNew = outs$cohortData, pixelGroupMapNew = outs$pixelGroupMap,
                         postFirePixelCohortData = postFirePixelCohortData,
                         burnedPixelCohortData, doAssertion = getOption("LandR.assertions", TRUE))
      
      sim$cohortData <- outs$cohortData
      sim$pixelGroupMap <- outs$pixelGroupMap
      sim$pixelGroupMap[] <- as.integer(sim$pixelGroupMap[])
    }
  }
  
  
  ## export objects
  sim$severityBMap <- severityBMap
  sim$severityData <- severityData
  sim$lastFireYear <- time(sim)
  
  ## TODO: Ceres: moved this to here to avoid re-killing/serotiny/repsoruting pixelGroups that burned in the previous year.
  ## update past pixelGroup number to match current ones.
  sim$treedFirePixelTableSinceLastDisp[, pixelGroup := as.integer(getValues(sim$pixelGroupMap))[pixelIndex]]
  # append previous year's
  treedFirePixelTableSinceLastDisp <- rbindlist(list(sim$treedFirePixelTableSinceLastDisp,
                                                     treedFirePixelTableSinceLastDisp))
  
  sim$treedFirePixelTableSinceLastDisp <- treedFirePixelTableSinceLastDisp
  
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  
  
  
  
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
  
  #  # input species ecoregion dynamics table
  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    sim$speciesEcoregion <- prepInputsSpeciesEcoregion(url = extractURL("speciesEcoregion"),
                                                       dPath = dPath, cacheTags = cacheTags)
  }
  
  
  return(invisible(sim))
}