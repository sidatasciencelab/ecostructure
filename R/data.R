#' Abundances of Himalayan bird species in local communities across elevations
#' A dataset of bird species abundances in local communities in the Himalayas
#' @format A matrix object with 38 sites (rows) and 
#' 304 species (columns) 
"him_birds_cn"

#' Morphological and functional metadata of Himalayan birds species
#' A dataset of morphological and functional traits of birds censused acrross
#' @format An ExpressionSet object with 38 sites (rows) and 
#' 304 species (columns) including phenoData which is site metadata and 
#' featureData which consists of feature metadata
"him_birds_md"

#' Local abundance data of Himalayan bird species with metadata
#' A dataset of bird abundance counts with site and feature metadata
#' @format An ExpressionSet object with 38 sites (rows) and 
#' 304 species (columns) including phenoData which is site metadata and 
#' featureData which consists of feature metadata
"him_grids_md"

#' An object of type phylo - which is a phylogenetic tree
#' @format An object of type phylo of the 304 bird species of the Himalayas
"phylo_tree"


#' The presence-absence data of birds from Australia
#' A dataset of bird presence absence data from Australia
#' @format A list comprising of two elements. One is 0/1 absence/presence matrix
#' of birds in the Australian continent - with sites along the rows and the 
#' bird species along the columns, with each entry being 0/1 based on
#' if the bird species is present in that site or not.  the other is a matrix
#' with equal number of rows as the presence -absence matrix with 2 columns 
#' representing the latitudes and longitudes.
"australia_birds"

#' The presence-absence data of birds from India/ S. Asia
#' A dataset of bird presence absence data from 
#' India/S. Asia
#' @format A list comprising of two elements. One is 0/1 absence/presence matrix
#' of birds in the Australian continent - with sites along the rows and the 
#' bird species along the columns, with each entry being 0/1 based on
#' if the bird species is present in that site or not.  the other is a matrix
#' with equal number of rows as the presence -absence matrix with 2 columns 
#' representing the latitudes and longitudes.
"indian_birds_pa"


#' The ecostructure_fit() model output on Australian birds
#' A list of output from K=6 GoM model run on Australian bird presence absence.
#' @format A list comprising of an omega and theta matrices obtained from
#' \code{ecostructure_fit()} on \code{australian_birds} data for K=6
"australia_model"


#' The ecostructure_fit() model on presence absence data from India/ S. Asia
#' A list of output from K=2 GoM model run on Indian bird presence absence.
#' @format A list comprising of an omega and theta matrices obtained from
#' \code{ecostructure_fit()} on \code{indian_birds_pa} data for K=2
"pres_ab_fit"


#' A dispersion field matrix example - used for creating a map
#' @format a dispersion field example created by \code{dispersion_fields_create}
#' in matrix format.
"dispersion_field_ex"
