#' `processPartialDist`
#'
#' @description
#' Processed partial disturbances (removal of cohorts) of pixels 
#'
#' @param cohortData
#' @param partialDistTable OPTIONAL need table or raster or both
#' @param partialDistRaster OPTIONAL need table or raster or both
#' @param inactivePixelIndex
#'
#' @return TODO
#'
#' @export
#' @importFrom data.table as.data.table
processPartialDist <- function(cohortData, partialDistTable, inactivePixelIndex) {
  
  #subset partialDistTable to disturbances of current year
  partialDistTable <- sim$partialDistTable[
    distYear == time(sim)]
  
  #create table of disturbed cohorts
  distCohorts <- sim$cohortData[
    partialDistTable,
    on = .(
      speciesCode,
      age >= ageMin,
      age <= ageMax
    ),
    nomatch = 0
  ]
  
  #create table with all cohorts on disturbed pixels
  colsToKeep <- c("speciesCode", "ecoregionGroup", "age", "B", "pixelGroup", "totalBiomass")
  postDistPixelCohortData <- sim$cohortData[
    distCohorts,
    on = .(pixelGroup),
    nomatch = 0
  ][, ..colsToKeep]
  
  #get disturbed pixel groups
  distPG <- unique(postDistPixelCohortData$pixelGroup)
  
  #get disturbed pixels
  distPixels <- which(values(sim$pixelGroupMap)[,1] %in% distPG)
  
  #build table of disturbed pixels and their associated pixel group.
  pixelLookup <- data.table(
    cellID     = distPixels,
    pixelGroup = values(sim$pixelGroupMap)[distPixels]
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
  
  #calcSiteShade requires a mortality column. Add it as NA
  #add check that only does this if the mortality column does not already exist
  postDistPixelTable[, mortality := NA]
  
  #site shade
  siteShade <- data.table(calcSiteShade(currentTime = round(time(sim)), postDistPixelTable,
                                        sim$speciesEcoregion, sim$minRelativeB))
  siteShade <- siteShade[, .(pixelGroup, siteShade)]
  postDistPixelTable <- siteShade[postDistPixelTable, on = "pixelGroup", nomatch = NA]
  postDistPixelTable[is.na(siteShade), siteShade := 0]
  rm(siteShade)
  
  #resprouting
  ##need a regular LandRCBM run to compare what all these tables are
  resproutingOutputs <- doResprouting(treedFirePixelTableSinceLastDisp = NA, #need to figure out my equivalent
                                      burnedPixelCohortData = postDistPixelTable,
                                      postFirePixelCohortData = postDistPixelCohortData,
                                      currentTime = time(sim), 
                                      species = sim$species,
                                      sufficientLight = sim$sufficientLight,
                                      calibrate = P(sim)$calibrate, #NULL
                                      postFireRegenSummary = sim$postFireRegenSummary) #NULL
  
  rePlanting <- plantCohorts(newPixelCohortData = distCohorts,
                             cohortData = sim$cohortData,
                             pixelGroupMap = sim$pixelGroupMap ,
                             initialB = 10,
                             currentTime = time(sim),
                             trackPlanting = FALSE
  )
  
  
  
  
  ## BUILD NEW COHORTS/PIXEL GROUPS
  
  # EXPORT OBJECTS
  ## new cohortData table
  ## new pixelGroupMap
  
}
