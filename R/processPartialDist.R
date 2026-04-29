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
processPartialDist <- function(cohortData, partialDistTable, pixelGroupMap, currentTime, partialDistLoc = NULL, ageMin = 50, ageMax = 500) {
  
  ## Subset partialDistTable to disturbances of current year and add age columns if not provided
  partialDistTable <- partialDistTable[distYear == currentTime]
  if (is.na(partialDistTable$ageMin)){
    partialDistTable$ageMin <- ageMin
  }
  if (is.na(partialDistTable$ageMax)){
    partialDistTable$ageMax <- ageMax
  }
  
  ## ID pixels that can be disturbed
  if (!is.null(partialDistLoc)) {
    if (inherits(partialDistLoc, "sf")) {
      partialDistLoc <- terra::vect(partialDistLoc)
      partialDistLoc <- terra::rasterize(
        partialDistLoc,
        pixelGroupMap,
        field = 1,     # all disturbed area gets value 1
        background = NA)}
    
    distValues <- terra::values(partialDistLoc, mat = FALSE)
    
    if (all(is.na(distValues))) { #if empty raster on which partialDistTable is to be applied
      coords <- crds(pixelGroupMap, df = TRUE)
      ext <- ext(partialDistLoc)
      
      distArea <- coords$x >= ext$xmin &
        coords$x <= ext$xmax &
        coords$y >= ext$ymin &
        coords$y <= ext$ymax
      
      distPixels <- which(distArea)
      
    } else { #if raster has values 1=disturbed 0=undisturbed
      distPixels <- which(distValues == 1)
    }
    
  } else {#partialDistLoc is NUlL and partialDistTable is applied over entire landscape
    distPixels <- seq_len(ncell(pixelGroupMap))
  }
  
  ## determine if all pixels in a pixelGroup is disturbed or now
  pixelGroups <- values(pixelGroupMap, mat = FALSE)
  pixelGroupsClean <- pixelGroups[!is.na(pixelGroups) & pixelGroups > 0]
  distPixelGroupsClean <- pixelGroups[distPixels]
  distPixelGroupsClean <- distPixelGroupsClean[!is.na(distPixelGroupsClean) & distPixelGroupsClean > 0]
  pixelGroupIDs <- sort(unique(pixelGroupsClean))
  totalTab <- tabulate(pixelGroupsClean, nbins = max(pixelGroupIDs))
  distTab <- tabulate(distPixelGroupsClean, nbins = max(pixelGroupIDs))
  totalPixelCount <- totalTab[pixelGroupIDs]
  distPixelCount <- distTab[pixelGroupIDs]
  
  pixelGroupDist <- data.table(
    pixelGroup = pixelGroupIDs,
    disturbedPixels = distPixelCount,
    totalPixels = totalPixelCount)
  
  pixelGroupDist[, disturbed := disturbedPixels > 0]
  pixelGroupDist[, partial := disturbedPixels > 0 & disturbedPixels < totalPixels]
  
  partialPG <- pixelGroupDist[partial == TRUE, pixelGroup]
  fullPG    <- pixelGroupDist[disturbed == TRUE & partial == FALSE, pixelGroup]
  
  ## create new pixel groups from groups that are not entirely disturbed (new group = disturbed, old group = undisturbed)
  maxPixelGroup <- max(cohortData$pixelGroup)
  # newPixelGroupmMap <- pixelGroups
  newPixelGroupmMap <- terra::values(pixelGroupMap, mat = FALSE, na.rm = FALSE)
  newCohorts <- list()
  newPGs <- integer()  #NEW
  
  for (pg in partialPG) {
    # create new pixel group
    maxPixelGroup <- maxPixelGroup + 1
    newPixelGroup <- maxPixelGroup
    
    newPGs <- c(newPGs, newPixelGroup) #NEW
    
    pgPixels <- which(pixelGroups == pg)
    disturbedSubset <- intersect(pgPixels, distPixels)
    
    # update pixelGroupMap with new pixelGroup in the disturbed area
    newPixelGroupmMap[disturbedSubset] <- newPixelGroup
    
    ## duplicate appropriate cohortData row for new disturbed group
    newRows <- copy(cohortData[pixelGroup == pg])
    newRows[, pixelGroup := newPixelGroup]
    
    newCohorts[[length(newCohorts) + 1]] <- newRows
  }
  
  ## update cohortData to add the new pixelGroup
  if (length(newCohorts) > 0)
    cohortData <- rbindlist(list(cohortData, rbindlist(newCohorts)))
  
  ## Update pixelGroupMap raster
  pixelGroupMap <- setValues(pixelGroupMap, newPixelGroupmMap)
  
  ## create table of disturbed cohorts
  # disturbedPixelGroups <- unique(newPixelGroupmMap[distPixels])
  disturbedPixelGroups <- unique(c(fullPG, newPGs)) #NEW
  distCohorts <- cohortData[pixelGroup %in% disturbedPixelGroups]
  distCohorts <- partialDistTable[
    distCohorts,
    on = .(speciesCode, ageMin <= age, ageMax >= age),
    nomatch = 0]
  
  ## clean up the columns to match cohortData
  distCohorts[, `:=`(
    age = ageMin,
    ageMin = NULL,
    ageMax = NULL)]
  colsToKeep <- c("speciesCode", "ecoregionGroup", "age", "B", "pixelGroup")
  distCohorts <- distCohorts[, ..colsToKeep]
  
  ## remove the distCohorts from sim$cohortData
  cohortData <- cohortData[!distCohorts, on = .(speciesCode, pixelGroup, age)]
  
  ## replanting
  rePlanting <- plantCohorts(
    newPixelCohortData = distCohorts,
    cohortData = cohortData,
    pixelGroupMap = pixelGroupMap,
    initialB = 10,
    currentTime = currentTime,
    trackPlanting = FALSE)
  
  ## recalc total biomass
  cohortData <- rePlanting[
    , totalBiomass := sum(B),
    by = pixelGroup]
  #ensure cohortData$pixelGroup is an integer
  cohortData[, pixelGroup := as.integer(pixelGroup)]
  ## return new cohortData and pixelGroupMap
  return(list(
    cohortData = cohortData,
    pixelGroupMap = pixelGroupMap
  ))
}