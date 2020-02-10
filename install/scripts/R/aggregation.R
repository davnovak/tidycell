## aggregation module

.get_label_vector <- function(file_idx, length_vector, junk = NULL) {
  labels            <- rep(file_idx, length_vector[file_idx])
  if (!is.null(junk)) {
    to_remove         <- unique(unlist(junk[[file_idx]]))
    labels[to_remove] <- -1
  }
  labels
}

## Function: convert vector of indices to list
.convert_to_list <- function(downsample.idcs.vector, # vector of sampling indices
                             lengths.vector) {       # vector of lengths of original samples
  
  list_idcs            <- cut(downsample.idcs.vector, breaks = c(1, cumsum(lengths.vector)), include.lowest = TRUE, labels = FALSE)
  downsample.idcs.list <- lapply(1:length(lengths.vector), function(idx) {
    which_idcs <- which(list_idcs == idx)
    if (length(which_idcs) == 0) return(integer(0))
    
    which_entries <- downsample.idcs.vector[which_idcs]
    
    if (idx > 1) which_entries <- which_entries - sum(lengths.vector[1:(idx - 1)]) # re-indexing
    which_entries
  })
  downsample.idcs.list # list of sampling indices, divided by sample of origin
}

## Function: using downsampling indices, collect data from .fcs files and create a concatenated downsampled matrix
.get_concat_matrix <- function(downsample_idcs,
                               file_paths) {
  
  concat_matrix <- lapply(1:length(downsample_idcs), function(idx) {
    file        <- file_paths[idx]
    which_rows  <- downsample_idcs[[idx]]
    read.FCS(file)@exprs[which_rows, ]
  }) %>% do.call(rbind, .)
}

## Function: create a concatenated expression matrix from multiple .fcs files
concatenate_fcs <- function(fcs.paths,
                            junk = NULL) {
  
  cat("Creating full concatenated expression matrix\n")
  efcs <- vector(mode = "list", length = length(fcs.paths))
  idcs <- vector(mode = "list", length = length(fcs.paths))
  for (idx in 1:length(fcs.paths)) {
    path <- fcs.paths[idx]
    file <- basename(path)
    ff   <- read.clean_FCS(path, junk = junk, id = file)
    
    if (idx == 1) cols <- colnames(ff)
    
    efcs[[idx]] <- ff@exprs
    if (idx == 1) {
      idcs[[idx]] <- c(1, nrow(ff@exprs))
    } else {
      previous_end <- idcs[[idx - 1]][2]
      idcs[[idx]]  <- c(previous_end + 1, previous_end + nrow(ff@exprs))
    }
    names(idcs)[idx] <- file
  }
  
  efcs <- do.call(rbind, efcs)
  
  colnames(efcs) <- cols
  
  cat("Concatenated expression matrix dimensions: ", nrow(efcs), " x ", ncol(efcs), "\n", sep = "")
  
  list(exprs       = efcs,
       sample_idcs = idcs)
}

## Function: create aggregate and balanced aggregate expression matrices from a collection of .fcs files
aggregate_data <- function(fcs.paths,
                           junk = NULL,
                           idcs.controls,
                           idcs.affected,
                           idcs.unused = NULL,
                           max_n_events = 5000000) {
  
  ## Get downsampling indices
  relative_idcs <- .reindex(idcs.controls, idcs.affected, idcs.unused) # reindex as if idcs.unused = {}
  
  if (!is.null(idcs.unused)) {
    paths <- fcs.paths[-idcs.unused]
    j     <- if (!is.null(junk)) { junk[-idcs.unused] } else { NULL }
  } else {
    paths <- fcs.paths
    j     <- if (!is.null(junk)) { junk } else { NULL }
  }
  rm(junk)
  
  len.all    <- sapply(paths, function(fcs) nrow(read.FCS(fcs)@exprs))
  labels.all <- unlist(lapply(1:length(paths), function(idx) .get_label_vector(file_idx = idx,
                                                                               length_vector = len.all,
                                                                               junk = j)))
  len.controls    <- len.all[relative_idcs$controls]
  labels.controls <- unlist(lapply(1:length(relative_idcs$controls), function(idx) .get_label_vector(file_idx = idx,
                                                                                                     length_vector = len.controls,
                                                                                                     junk = j[relative_idcs$controls])))
  len.affected    <- len.all[relative_idcs$affected]
  labels.affected <- unlist(lapply(1:length(idcs.affected), function(idx) .get_label_vector(file_idx = idx,
                                                                                            length_vector = len.affected,
                                                                                            junk = j[relative_idcs$affected])))
  N <- min(max_n_events, sum(labels.controls != -1), sum(labels.affected != -1))
  cat("Downsampling N =", N, "\n")
  
  cat("Computing downsampling indices\n")
  downsample_idcs.aggregate <- downsample(labels = labels.all, N = N, verbose = FALSE) %>% sort %>%
    .convert_to_list(., len.all)
  names(downsample_idcs.aggregate) <- basename(paths)
  
  downsample_idcs.controls <- downsample(labels = labels.controls, N = round(N/2), verbose = FALSE) %>% sort %>%
    .convert_to_list(., len.controls)
  names(downsample_idcs.controls) <- basename(paths[relative_idcs$controls])
  
  downsample_idcs.affected <- downsample(labels = labels.affected, N = round(N/2), verbose = FALSE) %>% sort %>%
    .convert_to_list(., len.affected)
  names(downsample_idcs.affected) <- basename(paths[relative_idcs$affected])
  
  ## Create aggregate matrix
  
  cat("Creating aggregate matrix\n")
  efcs.aggregate <- .get_concat_matrix(downsample_idcs.aggregate, paths)
  
  ## Create balanced aggregate matrix
  
  cat("Creating balanced aggregate matrix\n")
  efcs.controls  <- .get_concat_matrix(downsample_idcs.controls, fcs.paths[idcs.controls])
  efcs.affected  <- .get_concat_matrix(downsample_idcs.affected, fcs.paths[idcs.affected])
  efcs.balanced  <- rbind(efcs.controls, efcs.affected)
  
  list(aggregate                 = efcs.aggregate,
       balanced_aggregate        = efcs.balanced,
       controls_matrix           = efcs.controls,
       affected_matrix           = efcs.affected,
       downsample_idcs.aggregate = downsample_idcs.aggregate,
       downsample_idcs.controls  = downsample_idcs.controls,
       downsample_idcs.affected  = downsample_idcs.affected)
}
