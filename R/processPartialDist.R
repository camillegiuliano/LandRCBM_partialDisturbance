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
  postDistPixelCohortData <- sim$cohortData[
    distCohorts,
    on = .(pixelGroup),
    nomatch = 0
  ]
  
  # id disturbed pixels
  distPixelMap <- classify(
    sim$pixelGroupMap,
    rcl = cbind(unique(distCohorts$pixelGroup), 1),
    others = NA
  )
  
  #get IDs of disturbed pixels"
  distPixels <- which(values(sim$pixelGroupMap)[,1] %in% unique(distCohorts$pixelGroup))
  
  
  
  
  
  
}
















