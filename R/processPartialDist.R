#' `processPartialDist`
#'
#' @description
#' Processed partial disturbances (removal of cohorts) of pixels 
#'
#' @param cohortData
#' @param PartialDistTable
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
  
  #site shade
  siteShade <- data.table(calcSiteShade(currentTime = round(time(sim)), postDistPixelTable,
                                        sim$speciesEcoregion, sim$minRelativeB))
  ##MORTALITY COLUMN
  
  #resprouting
  
}
