## pre-processing module

.is_single <- function(ff, plot = FALSE, nMAD = 4) {
  
  require(flowCore)
  fsc_a  <- exprs(ff)[,"FSC-A"]
  fsc_h  <- exprs(ff)[,"FSC-H"]
  ratios <- fsc_h / fsc_a
  
  tr        <- median(ratios) - nMAD * mad(ratios)
  selection <- ratios > tr
}

## Function: set lower and upper limit in each channel; anything past these values is collapsed to the extreme (eg. for every x > maximum: x <- maximum).
## If out.idcs is not NULL, assign the indices of these values which exceeded the minimum/maximum (for later removal of events which are
## extreme in at least one non-scatter channel).
.set_signal_lim <- function(ff, which_channels, y.min, y.max, out.idcs = NULL, marker_names = NULL) {
  e <- ff@exprs[, which_channels]
  if (is.null(marker_names)) marker_names <- which_channels
  
  channels            <- colnames(e)
  scatter_channels    <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")
  channels.no_scatter <- channels[which(!channels %in% scatter_channels)]
  
  idcs <- list()
  if (!is.null(y.min)) {
    if (!is.null(out.idcs)) {
      extremes  <- apply(e[, channels.no_scatter], 1, function(y) any(y < y.min))
      idcs[[1]] <- extremes
    }
    e[e < y.min] <- y.min 
  }
  if (!is.null(y.max)) {
    if (!is.null(out.idcs)) {
      extremes  <- apply(e[, channels.no_scatter], 1, function(y) any(y > y.max))
      idcs[[2]] <- extremes
    }
    e[e > y.max] <- y.max
  }
  ff@exprs[, which_channels] <- e
  
  
  
  if (!is.null(out.idcs)) eval.parent(substitute(out.idcs <- unique(unlist(idcs))))
  ff
}

## Function: like .set_signal_lim, but allows writing a file identifying the channels in which extreme values occur
.set_signal_lim_and_report <- function(ff, which_channels, y.min, y.max, out.idcs = NULL, marker_names = NULL, output_file = NULL) {
  ## TO DO: rewrite to use a stream for writing the output file
  e <- ff@exprs[, which_channels]
  if (is.null(marker_names)) marker_names <- which_channels
  idcs <- list()
  if (!is.null(y.min)) {
    if (!is.null(out.idcs)) {
      violations <- apply(e, 1, function(y) which(y < y.min))
      extremes   <- which(lapply(violations, length) > 0)
      idcs[[1]]  <- extremes
      
      if (!is.null(output_file)) cat("Low extremes...\n", file = output_file, append = TRUE)
      
      for (i in violations) {
        if (length(i) > 0) {
          text <- paste0(paste0(marker_names[i], collapse = ", "), "\n")
          if (!is.null(output_file)) cat(text, file = output_file, append = TRUE)
        }
      }
    }
    e[e < y.min] <- y.min 
  }
  if (!is.null(output_file)) cat("", file = output_file, append = TRUE)
  
  if (!is.null(y.max)) {
    if (!is.null(out.idcs)) {
      violations <- apply(e, 1, function(y) which(y > y.max))
      extremes   <- which(lapply(violations, length) > 0)
      idcs[[2]]  <- extremes
      
      if (!is.null(output_file))  cat("High extremes...", file = output_file, append = TRUE)
      
      for (i in violations) {
        if (length(i) > 0) {
          text <- paste0(paste0(marker_names[i], collapse = ", "), "\n")
          if (!is.null(output_file)) cat(text, file = output_file, append = TRUE)
        }
      }
    }
    e[e > y.max] <- y.max
  }
  ff@exprs[, which_channels] <- e
  
  if (!is.null(out.idcs)) eval.parent(substitute(out.idcs <- unique(unlist(idcs))))
  
  ff
}

## Function: write pre-processed .fcs files (clipping extreme values, compensation, transformation) and identify
## junk events for later removal
preprocess_fcs <- function(fcs.paths,
                           to_remove_immediately = NULL, # discouraged
                           known_junk = NULL,
                           scatter_idcs = NULL,
                           compensation_matrix = NULL,
                           adjust_signal_extremes = NULL,
                           flag_doublets = FALSE,
                           flag_debris = FALSE,
                           low_quality_metric = NULL, # using flowAI
                           transformation = NULL,
                           dead_cells.channel = NULL,
                           dead_cells.threshold = 1.7,
                           dead_cells.threshold.tolerance = 0.3,
                           foldername = "fcs_files_preprocessed",
                           out.events_to_remove = NULL,
                           make_plots = FALSE,
                           visual.path = NULL) {
  
  if (!is.null(visual.path)) {
    d <- ceiling(sqrt(length(fcs.paths) * (2 + (!is.null(dead_cells.channel)))))
    png(visual.path,
        width  = d * 600,
        height = d * 600)
    layout(matrix(1:(d * d), ncol = d, byrow = TRUE))
  }
  
  if (!file.exists(foldername)) dir.create(foldername)
  
  cat("Initial pre-processing...\n")
  
  junk_idcs <- vector(mode = "list", length = length(fcs.paths)); names(junk_idcs) <- basename(fcs.paths)
  for (i in 1:length(fcs.paths)) {
    file <- fcs.paths[i]
    cat("(", i, "/", length(fcs.paths), ") ", file, "\n", sep = "")
    ff   <- read.FCS(file, truncate_max_rage = FALSE)
    
    junk_extremes <- NULL
    
    if (!is.null(adjust_signal_extremes)) {
      which.channels <- colnames(ff@exprs)
      which.channels <- which.channels[!which.channels %in% adjust_signal_extremes$exclude_channels]
      junk_extreme   <- 0
      ff             <- .set_signal_lim(ff = ff,
                                        which_channels = which.channels,
                                        y.min = adjust_signal_extremes$min,
                                        y.max = adjust_signal_extremes$max,
                                        out.idcs = junk_extreme)
    }
    
    if (!is.null(compensation_matrix)) ff <- compensate(ff, compensation_matrix)
    
    if (!is.null(transformation)) {
      if (is.character(transformation)) {
        if (transformation == "estimateLogicle") {
          message('Applying estimageLogicle transformation')
          channels.transform <- colnames(ff)
          if (!is.null(scatter_idcs)) channels.transform <- channels.transform[-scatter_idcs]
          channels.transform <- channels.transform[!channels.transform %in% c("Time", "TIME")]
          lgcl <- estimateLogicle(ff, channels.transform)
          ff   <- flowCore::transform(ff, lgcl)
        }
      } else {
        message('Applying asinh transformation')
        ff <- flowCore::transform(ff, transformation)
      }
      if (!is.null(scatter_idcs)) ff@exprs[, scatter_idcs] <- scale(ff@exprs[, scatter_idcs])
    }
    
    junk_debris <- NULL
    if (make_plots) plotDens(ff, c("FSC-A", "SSC-A"), main = basename(file))
    if (flag_debris) {
      debris_cuts <- deGate(ff, "FSC-A", all.cuts = TRUE, upper = FALSE, verbose = FALSE)
      debris_cut  <- debris_cuts[which.min(abs(debris_cuts - 30000))]
      junk_debris <- which(ff@exprs[, "FSC-A"] <= debris_cut)
      
      abline(v = debris_cuts, col = "grey")
      abline(v = debris_cut, col = "red", lwd = 2)
    }
    
    junk_doublets <- NULL
    if (make_plots) plotDens(ff, c("FSC-A", "FSC-H"), main = basename(file))
    if (flag_doublets) {
      junk_doublets <- which(!(.is_single(ff, nMAD = 5)))
      points(ff@exprs[junk_doublets, c("FSC-A", "FSC-H")], col = "red")
    }
    
    junk_dead <- NULL
    if (!is.null(dead_cells.channel)) {
      plotDens(ff, c(dead_cells.channel, "FSC-A"), main = basename(file))
      dead_cuts <- deGate(ff, dead_cells.channel, all.cuts = TRUE, upper = FALSE, verbose = FALSE)
      dead_cut  <- dead_cuts[which.min(abs(dead_cuts - dead_cells.threshold))]
      if (abs(dead_cut - dead_cells.threshold) > dead_cells.threshold.tolerance) dead_cut <- dead_cells.threshold
      abline(v = dead_cut, col = "red", lwd = 2)
      junk_dead <- which(ff@exprs[, dead_cells.channel] > dead_cut)
    }
    
    if (!is.null(out.events_to_remove)) {
      junk_idcs[[i]]          <- list()
      junk_idcs[[i]]$debris   <- junk_debris
      junk_idcs[[i]]$doublets <- junk_doublets
      junk_idcs[[i]]$dead     <- junk_dead
      junk_idcs[[i]]$extreme  <- junk_extremes
      junk_idcs[[i]]$known    <- known_junk[[i]]
    }
    
    if (!is.null(to_remove_immediately)) {
      if (length(to_remove_immediately[[i]]) > 0) ff@exprs <- ff@exprs[-to_remove_immediately[[i]], ]
    }
    
    write_path <- file.path(foldername, basename(file))
    write.FCS(make_valid_fcs(ff@exprs, ff), write_path)
  }
  
  if (!is.null(visual.path)) dev.off()
  
  if (!is.null(low_quality_metric)) {
    cat("Quality control using flowAI...\n")
    for (i in 1:length(fcs.paths)) {
      file              <- file.path(foldername, basename(fcs.paths[i]))
      cat("(", i, "/", length(fcs.paths), ") ", file, "\n", sep = "")
      ff                <- read.FCS(file)
      anomalous         <- flow_auto_qc(ff, output = 3, remove_from = low_quality_metric)
      gc()
      junk_idcs[[i]]$qc <- unlist(anomalous)
    }
  }
  
  if (!is.null(out.events_to_remove)) eval.parent(substitute(out.events_to_remove <- junk_idcs))
}

