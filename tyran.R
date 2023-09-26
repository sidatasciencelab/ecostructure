any.allzero.cols(pres_abs_global_sparse)
hbw_birds_taxonomy <- as.matrix(read.csv("data/Handbook of the Birds of the World and BirdLife International Digital Checklist of the Birds of the World_Version_6b.csv"))
hbw_birds_taxonomy <- hbw_birds_taxonomy[-c(1,2),]
colnames(hbw_birds_taxonomy) <- hbw_birds_taxonomy[1,]
hbw_birds_taxonomy <- hbw_birds_taxonomy[-1,]
hbw_birds_taxonomy <- as.data.frame(hbw_birds_taxonomy) %>% filter(`BirdLife taxonomic treatment:\nR = recognised as a species;\nNR = not recognised as a species` == "R")
hbw_tyrannidea <- hbw_birds_taxonomy %>% filter(`Family name` == "Tyrannidae")
hbw_tyrannidae <- hbw_tyrannidea
save(hbw_tyrannidea, file ="data/hbw_tyrannidae.rda")




dsp_to_matrix <- function(dispersion.field) {
  cat("\n Converting dispersion fields to matrix rows \n")
  pb <- txtProgressBar()
  temp_data <- dispersion.field[[1]]
  temp_data[is.na(temp_data)] <- 0
  map_data <- data.frame(y = rep(NA, (dim(temp_data[[1]])[1] * dim(temp_data[[1]])[2])), x = rep(NA, length(temp_data)))
  for (l in 1:length(dispersion.field)) {
    temp_data <- dispersion.field[[l]]
    temp_data[is.na(temp_data)] <- 0
    temp_data_vec <- matrix(as.matrix(temp_data),
                            ncol = dim(temp_data)[1] * dim(temp_data)[2], byrow = TRUE)
    map_data[l, ] <- temp_data_vec
    setTxtProgressBar(pb, l / length(dispersion.field))
  }
  rownames(map_data) <- names(dispersion.field)
  return(map_data)
}






create <- function(elems) {
  return(as.data.table(elems))
}

append <- function(dt, elems) {
  n <- attr(dt, "rowcount")
  if (is.null(n)) {
    n <- nrow(dt)
  }
  if (n == nrow(dt)) {
    tmp <- elems[1]
    tmp[[1]] <- rep(NA, n)
    dt <- rbindlist(list(dt, tmp), fill = TRUE, use.names = TRUE)
    setattr(dt, "rowcount", n)
  }
  pos <- as.integer(match(names(elems), colnames(dt)))
  for (j in seq_along(pos))
  {
    set(dt, i = as.integer(n + 1), pos[[j]], elems[[j]])
  }
  setattr(dt, "rowcount", n + 1)
  return(dt)
}

access <- function(elems) {
  n <- attr(elems, "rowcount")
  return(as.data.table(elems[1:n, ]))
}


dsp_from_scratch$raster