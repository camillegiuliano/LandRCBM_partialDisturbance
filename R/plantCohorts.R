#' `plantCohorts`
#'
#' @description
#' Processed partial disturbances (removal of cohorts) of pixels 
#'
#' @param newPixelCohortData
#' @param cohortData
#' @param pixelGroupMap 
#' @param initialB default is 10
#' @param currentTime
#' @param successionTimestep
#' @param trackPlanting default is FALSE
#'
#' @return TODO
#'
#' @export
#' @importFrom data.table as.data.table
plantCohorts <- function (newPixelCohortData, cohortData, pixelGroupMap, initialB = 10, 
                          currentTime, successionTimestep, trackPlanting = FALSE) 
{
  namesNCD <- names(newPixelCohortData)
  if (!isTRUE("pixelGroup" %in% namesNCD)) {
    if (isTRUE("pixelIndex" %in% namesNCD)) {
      newPixelCohortData[, `:=`(pixelGroup, as.vector(pixelGroupMap[])[pixelIndex])]
    }
    else {
      stop("newPixelCohortData must have either pixelIndex or pixelGroup")
    }
  }
  if (!is.null(newPixelCohortData[["pixelIndex"]])) {
    set(newPixelCohortData, NULL, "pixelIndex", NULL)
  }
  duplicates <- duplicated(newPixelCohortData[, .(speciesCode, 
                                                  ecoregionGroup, pixelGroup, age)])
  newCohortData <- newPixelCohortData[!duplicates]
  newCohortData[, `:=`(age, 2)]
  if (isTRUE(is.na(initialB)) || is.null(initialB)) {
    initialB <- formals(plantNewCohorts)$initialB
  }
  newCohortData[, `:=`(B, initialB)]
  newCohortData <- newCohortData[, .(pixelGroup, ecoregionGroup, 
                                     speciesCode, age, B, mortality = 0L, aNPPAct = 0L)]
  if (getOption("LandR.assertions")) {
    if (isTRUE(NROW(unique(newCohortData, by = c("pixelGroup", 
                                                 "age", "speciesCode"))) != NROW(newCohortData))) {
      stop("Duplicated new cohorts in a pixelGroup. Please debug LandR:::.plantNewCohorts")
    }
  }
  if (trackPlanting) {
    newCohortData[, `:=`(planted, TRUE)]
  }
  cohortData <- rbindlist(list(cohortData, newCohortData), 
                          fill = TRUE, use.names = TRUE)
  return(cohortData)
}
