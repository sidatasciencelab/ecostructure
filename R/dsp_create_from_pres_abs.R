#' @title Create dispersion fields and a presence absence matrix from geodatabase 
#'        file or a shape file containing all species.
#'
#' @description Reads a single GIS data file and generates a presence-matrix
#' and associated dispersion fields in a user-specified region.
#'
#' @param gdb_object The geodatabase or a single shapefile containing all species to be read 
#'                    in ahead of time and provided as an sf object. 
#' @param raster_latlim The latitudinal range of the desired dispersion fields and presence-absence matrix.
#' @param raster_longlim The longitudinal range.
#' @param raster_resolution The resolution of the raster in degrees.
#' @param precision The precision of the initial estimation provided to fasterize. 
#'                          fasterize tends to undercount small species ranges, so is required 
#'                          to keep at a low value in order detect species presences.
#'                          Most users will not need to change, but increaseing this value will 
#'                          speed up the function but with potential for error. 
#' @param thresh The number of co-occuring species required to generate a dispersion field for that map cell. 
#' @param drop Boolean. If T sites with 0 overlapping species will be dropped from the final presence-absence matrix.
#'             Defaults to false.  
#' @param species_feature The name of the feature (column) that contains the species names 
#'                             in the GIS data source.
#' @param species_names If you would like to only include a subset of the species is the GIS data source
#'                      a vector of species names can be provided to filter the GIS data before processing
#'                      into a presence absence matrix and dispersion fields. 
#'
#' @return Returns 2 lists and 1 matrix. In each list, each element is a
#' dispersion field corresponding to the species present in the 
#' specific site to which it corresponds. The first list "matrix" 
#' contains each dispersion field as a matrix. The second list "raster" returns 
#' each dispersion field as a raster. The matrix is a presence absence matrix with 
#' species as columns and each individual map cell the rows. The names of the elements 
#' of the lists, and the rownames for the presence-absence matrix, are assigned as 
#' the lat_long combos for those locations.
#'
#' @keywords counts data, local site-species data, dispersion field
#'
#' @import sf
#' @importFrom fasterize fasterize
#' @importFrom raster raster setValues stack aggregate as.matrix
#' @importFrom dplyr filter
#' @export


dsp_create_from_pres_abs <- function(pres_ab,
                                     r_sfs_ag,
                                     raster_resolution = 5,
                                     thresh = 4,
                                     raster_latlim = c(-90, 90),
                                     raster_longlim = c(-180, 180),
                                     species_feature = "SCINAME",
                                     precision = 0.25,
                                     species_names = NULL,
                                     drop = FALSE) {
  
  r <- raster::raster(resolution = raster_resolution, xmn = raster_longlim[1], 
                      xmx = raster_longlim[2], ymn = raster_latlim[1], ymx = raster_latlim[2])
  
  r_prim <- raster::raster(resolution = precision, xmn = raster_longlim[1], 
                           xmx = raster_longlim[2], ymn = raster_latlim[1], 
                           ymx = raster_latlim[2])
  
  idx <- list()
  for (i in 1:raster::ncell(r)) {
    idx[[i]] <- as.vector(raster::xyFromCell(r, i, spatial = FALSE))
  }
  
  overlps_sp <- apply(pres_ab, 2, function(x) which(x == 1))
  overlps_cell <- apply(pres_ab, 1, function(x) which(x == 1))

  cells_with_overlap <- which(unlist(lapply(overlps_cell, any)))
  focal_LL <- idx[cells_with_overlap]
  cells_with_asblg <- which(unlist(lapply(overlps_cell, function(x) length(x) > thresh)))
  focal_asblg_LL <- idx[cells_with_asblg]

  dispersion.field.raster <- list()
  dispersion.field.matrix <- list()

  lat_long_names <- c()

  cat("\n Generating disersion fields \n")
  pb <- txtProgressBar()
  for (i in 1:length(cells_with_asblg)) {
    dsp_fld <- sum(stack(r_sfs_ag[overlps_cell[cells_with_asblg[i]][[1]]]))
    setTxtProgressBar(pb, i / length(cells_with_asblg))
    dsp_fld[dsp_fld == 0] <- NA
    dispersion.field.raster[[i]] <- dsp_fld
    dispersion.field.matrix[[i]] <- raster::as.matrix(dsp_fld)
    lat_long_names[i] <- paste(focal_asblg_LL[[i]][2],
      focal_asblg_LL[[i]][1],
      sep = "_"
    )
  }

  names(dispersion.field.matrix) <- lat_long_names
  names(dispersion.field.raster) <- lat_long_names

  if (drop) {
    pres_ab <- pres_ab[rowSums(pres_ab) > 0, ]
  }

  dispersion.field <- list(
    matrix = dispersion.field.matrix,
    raster = dispersion.field.raster,
    pres_ab = pres_ab
  )

  return(dispersion.field)
}
