#' `processPartialDist`
#'
#' @description
#' Processed partial disturbances (removal of cohorts) of pixels 
#'
#' @param cohortData
#' @param partialDistTable OPTIONAL need table or raster or both
#' @param partialDistRaster OPTIONAL need table or raster or both
#' @param pixelGroupMap
#'
#' @return TODO
#'
#' @export
#' @importFrom data.table as.data.table
processPartialDist <- function(cohortData, partialDistTable, pixelGroupMap, currentTime, partialDistLoc = NULL) {
  
  #subset partialDistTable to disturbances of current year
  partialDistTable <- partialDistTable[
    distYear == currentTime]
  
  #create table of disturbed cohorts
  distCohorts <- partialDistTable[
    cohortData,
    on = .(
      speciesCode,
      ageMin <= age,
      ageMax >= age),
    nomatch = 0]
  
  #clean up the columns to match cohortData
  distCohorts[
    , `:=`(
      age = ageMin,   
      ageMin = NULL, 
      ageMax = NULL)]
  colsToKeep <- c("speciesCode", "ecoregionGroup", "age", "B", "pixelGroup")
  distCohorts <- distCohorts[, ..colsToKeep]
  
  
  #create table with all cohorts on disturbed pixels
  colsToKeep <- c("speciesCode", "ecoregionGroup", "age", "B", "pixelGroup")
  postDistPixelCohortData <- cohortData[
    distCohorts,
    on = .(pixelGroup),
    nomatch = 0
  ][, ..colsToKeep]
  
  #get disturbed pixel groups
  distPG <- unique(postDistPixelCohortData$pixelGroup)
  
  #get disturbed pixels
  distPixels <- which(values(pixelGroupMap)[,1] %in% distPG)
  
  #build table of disturbed pixels and their associated pixel group.
  pixelLookup <- data.table(
    cellID     = distPixels,
    pixelGroup = values(pixelGroupMap)[distPixels]
  )
  
  #pixel-specific table with all cohort information
  postDistPixelTable <- postDistPixelCohortData[
    pixelLookup,
    on = .(pixelGroup),
    nomatch = 0,
    allow.cartesian = TRUE
  ]
  
  #remove disturbed cohorts and recalc totalBiomass
  spRemove <- unique(partialDistTable$speciesCode)
  postDistPixelTable <- postDistPixelTable[
    !speciesCode %in% spRemove
  ][
    , totalBiomass := sum(B),
    by = cellID
  ]
  
  postDistPixelCohortData <- postDistPixelCohortData[
    !speciesCode %in% spRemove
  ]
  
  
  #calcSiteShade requires a mortality column. Add it as NA
  #add check that only does this if the mortality column does not already exist
  postDistPixelTable[, mortality := NA]
  
  # #site shade
  # siteShade <- data.table(calcSiteShade(currentTime = round(time(sim)), postDistPixelTable,
  #                                       sim$speciesEcoregion, sim$minRelativeB))
  # siteShade <- siteShade[, .(pixelGroup, siteShade)]
  # postDistPixelTable <- siteShade[postDistPixelTable, on = "pixelGroup", nomatch = NA]
  # postDistPixelTable[is.na(siteShade), siteShade := 0]
  # rm(siteShade)
  # 
  # #resprouting
  # ##need a regular LandRCBM run to compare what all these tables are
  # resproutingOutputs <- doResprouting(treedFirePixelTableSinceLastDisp = NA, #need to figure out my equivalent
  #                                     burnedPixelCohortData = postDistPixelTable,
  #                                     postFirePixelCohortData = postDistPixelCohortData,
  #                                     currentTime = time(sim), 
  #                                     species = sim$species,
  #                                     sufficientLight = sim$sufficientLight,
  #                                     calibrate = P(sim)$calibrate, #NULL
  #                                     postFireRegenSummary = sim$postFireRegenSummary) #NULL
  
  
  #remove the distCohorts from sim$cohortData
  cohortData <- cohortData[!distCohorts, on =.(speciesCode, pixelGroup, age)]
  
  rePlanting <- plantCohorts(newPixelCohortData = distCohorts,
                             cohortData = cohortData,
                             pixelGroupMap = pixelGroupMap ,
                             initialB = 10,
                             currentTime = currentTime,
                             trackPlanting = FALSE
  )
  #FIGURE OUT how to make work with out the provenance table. What is the provenance table??
  # postDistCohortTable <- updateCohortData(newPixelCohortData = rePlanting,
  #                                                    cohortData = cohortData,
  #                                                    pixelGroupMap = sim$pixelGroupMap,
  #                                                    currentTime = time(sim),
  #                                                    speciesEcoregion = sim$speciesEcoregion,
  #                                                    initialB = 10,
  #                                                    cohortDefinitionCols = LandR::cohortDefinitionCols(),
  #                                                    verbose = getOption("LandR.verbose", TRUE),
  #                                                    doAssertion = getOption("LandR.assertions", TRUE)
  # )
  
  #recalc totalBiomass here
  cohortData <- rePlanting[
    , totalBiomass := sum(B),
    by = pixelGroup
  ]
  
  ## BUILD NEW COHORTS/PIXEL GROUPS
  
  # EXPORT OBJECTS
  ## new cohortData table
  ## new pixelGroupMap
  return(cohortData)
}