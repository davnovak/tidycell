## population thresholding module

population_thresholds <- function(fcs.paths,
                                  junk,
                                  channels_of_interest = NULL,
                                  markers_of_interest = NULL,
                                  default_thresholds = NULL,
                                  tolerance = NULL,
                                  visual.folder = NULL) {
  
  if (!is.null(visual.folder)) {
    if (!file.exists(visual.folder)) dir.create(visual.folder)
    d <- ceiling(sqrt(length(channels_of_interest)))
  }
  
  thresholds <- vector(mode = "list", length = length(fcs.paths)); names(thresholds) <- fcs.paths
  for (i in 1:length(fcs.paths)) {
    file <- fcs.paths[i]; cat("(", i, "/", length(fcs.paths), ") ", file, "\n", sep = "")
    ff   <- read.clean_FCS(fcs = file, junk = junk, id = basename(file))
    
    threshold <- lapply(channels_of_interest, function(channel) {
      cuts    <- suppressWarnings(deGate(ff, channel, all.cuts = TRUE, upper = FALSE, verbose = FALSE))
      default <- if (length(default_thresholds) == 1) { default_thresholds[[1]] } else { default_thresholds[[channel]] }
      
      if (is.na(cuts) || length(cuts) == 0) {
        cuts <- default
      } else {
        cuts <- cuts[which.min(abs(cuts - default))][1]
        if (abs(cuts - default) > tolerance) cuts <- default
      }
      
      cuts
    })
    
    threshold        <- unlist(threshold)
    names(threshold) <- markers_of_interest
    thresholds[[i]]  <- threshold
    
    if (!is.null(visual.folder)) {
      png(filename = file.path(visual.folder, paste0(strsplit(basename(file), "[.]")[[1]][1], "_thresholds.png")),
          width    = d * 500,
          height   = d * 500)
      layout(matrix(1:(d*d), ncol = d, byrow = TRUE))
      for (j in 1:length(channels_of_interest)) {
        channel <- channels_of_interest[j]
        marker  <- markers_of_interest[j]
        plotDens(ff, c(channel, "FSC-A"), main = "Threshold")
        abline(v = threshold[marker], lwd = 3, lty = 2, col = "red")
      }
      dev.off()
    }
  }
  
  thresholds
}