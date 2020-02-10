## input filtering module

## Function: read .fcs file and exclude junk events from its expression matrix
read.clean_FCS <- function(fcs,
                           junk,
                           id = NULL) {
  
  ff <- read.FCS(fcs, truncate_max_range = FALSE)
  if (is.null(junk)) return(ff)
  to_remove <- unique(unlist(junk[[id]]))
  if (length(to_remove) > 0) ff@exprs <- ff@exprs[-to_remove, ]
  
  ff
}

## Function: equivalent of read.clean_FCS for vectors of labels per event
clean_vector <- function(vec,
                         junk,
                         id = NULL) {
  
  if (is.null(junk)) return(vec)
  to_remove <- unique(unlist(junk[[idx]]))
  vec       <- vec[-to_remove]
}