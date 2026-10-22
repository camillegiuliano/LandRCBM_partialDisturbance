## processPartialDist(): removal of cut cohorts, and replanting
library(data.table)
library(terra)
library(LandR)
for (f in c("processPartialDist.R", "plantCohorts.R")) source(file.path("..", "..", "R", f))

test_that("a species with two cut cohorts in a pixelGroup is replanted once", {
  ## one pixelGroup, cut everywhere; fir at 50 and 63 (both cuttable), spruce at 20 (too young)
  pixelGroupMap <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, vals = 1L)
  cohortData <- data.table(pixelGroup = 1L, ecoregionGroup = factor("1_01"),
                           speciesCode = factor(c("Abie_las", "Abie_las", "Pice_eng")),
                           age = c(50L, 63L, 20L), B = c(50L, 4356L, 300L))
  partialDistTable <- data.table(distYear = 2075, speciesCode = c("Abie_las", "Pice_eng"))
  out <- processPartialDist(cohortData, partialDistTable, pixelGroupMap, replanting = TRUE,
                            currentTime = 2075)
  cd <- out$cohortData
  expect_false(any(cd$speciesCode == "Abie_las" & cd$age >= 50))  # both fir cohorts removed
  expect_equal(nrow(cd[speciesCode == "Abie_las"]), 1)            # replanted once
  expect_equal(cd[speciesCode == "Pice_eng", age], 20)              # too young: untouched
})
