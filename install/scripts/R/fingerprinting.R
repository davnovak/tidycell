## fingerprinting module

.sgn <- function(x) ifelse(x >= 0, 1, -1) # signum

aberrance <- function(aa_s3,
                      fcs.paths,
                      junk             = NULL,
                      idcs.controls,
                      idcs.affected,
                      idcs.unused      = NULL,
                      channels_of_interest,
                      markers_of_interest,
                      method,
                      method.settings,
                      tidycell_path    = "",
                      n_bins           = 32,
                      output_directory = NULL) {
  
  ## Reindex as if idcs.unused = {}
  idcs.controls.original <- idcs.controls
  idcs.affected.original <- idcs.affected
  fcs.paths.original     <- fcs.paths
  
  if (!is.null(idcs.unused)) {
    fcs.paths     <- fcs.paths[-idcs.unused]
    relative_idcs <- .reindex(idcs.controls, idcs.affected, idcs.unused)
    idcs.controls <- relative_idcs$controls
    idcs.affected <- relative_idcs$affected
  }
  
  cat("Computing fingerprints of samples. The following files are included in this iteration:", paste0(basename(fcs.paths), collapse = "\n"), sep = "\n")
  
  if (!is.null(idcs.unused))
    cat("The following files are excluded from this iteration:", paste0(basename(fcs.paths[-idcs.unused]), collapse = "\n"), sep = "\n")
  
  bins <- vector(mode = "list", length = length(fcs.paths)); names(bins) <- basename(fcs.paths)
  bin_overlap <- FALSE
  
  if (method == "hyperbinning") { ###############################
    
    cat("Method: hyperbinning\n")
    
    if (log2(n_bins)%%1 != 0) stop("For hyperbinning, number of bins must be a power of 2")
    
    cat("Number of bins:", n_bins, "\n")
    
    model <- HBmodel(flowFrame(aa_s3$agg$balanced_aggregate[, which(colnames(aa_s3$agg$balanced_aggregate) %in% channels_of_interest)]),
                     depth = log2(n_bins),
                     channels = channels_of_interest)
    
    for (i in 1:length(fcs.paths)) {
      file <- fcs.paths[i]
      ff   <- read.clean_FCS(file, junk = junk, id = basename(file))
      
      # If channels mismatched, match by markers
      ff_channels <- colnames(ff@exprs)
      if (!all(channels_of_interest %in% ff_channels)) {
        ff_markers <- markernames(ff)
        
        if (!all(markers_of_interest %in% ff_markers)) stop(paste0("Cannot find channels/markers of interest in ", file))
        
        channel_idcs       <- which(markers_of_interest %in% ff_markers)
        ff@exprs           <- ff@exprs[, channel_idcs]
        colnames(ff@exprs) <- channels_of_interest
      }
      
      classification       <- HBclassify(model, ff)
      
      bins[[i]]$assignment <- classification$visualize
      fp                   <- classification$finger
      bins[[i]]$fp.abs     <- fp           # non-scaled fingerprint
      bins[[i]]$fp.rel     <- fp / sum(fp) # proportionate
    }
    names(bins) <- basename(fcs.paths)
    
    aggregate_matrix_binning          <- HBclassify(model, flowFrame(aa_s3$agg$aggregate))$visualize
    balanced_aggregate_matrix_binning <- HBclassify(model, flowFrame(aa_s3$agg$balanced_aggregate))$visualize
    
  } else if (method == "FlowSOM") { ###############################
    
    cat("Method: FlowSOM\n")
    
    cat("Number of bins (SOM metaclusters):", n_bins, "\n")
    
    scale  <- if (!is.null(method.settings$scale))  { method.settings$scale  } else { FALSE }
    grid.x <- if (!is.null(method.settings$grid.x)) { method.settings$grid.x } else { 15 }
    grid.y <- if (!is.null(method.settings$grid.y)) { method.settings$grid.y } else { 15 }
    
    som.input <- ReadInput(flowFrame(aa_s3$agg$balanced_aggregate), scale = scale, silent = TRUE)
    som.grid  <- BuildSOM(som.input, colsToUse = which(colnames(aa_s3$agg$balanced_aggregate) %in% channels_of_interest), xdim = method.settings$grid.x, ydim = method.settings$grid.y, silent = TRUE)
    som.grid  <- BuildMST(som.grid, silent = TRUE)
    som.grid$prettyColnames[names(som.grid$prettyColnames) %in% channels_of_interest] <- markers_of_interest
    
    som.metaClustering <- metaClustering_consensus(som.grid$map$codes, k = n_bins)
    
    for (i in 1:length(fcs.paths)) {
      file <- fcs.paths[i]
      ff   <- read.clean_FCS(file, junk = junk, id = basename(file))
      
      # If channels mismatched, match by markers
      ff_channels <- colnames(ff@exprs)
      if (!all(channels_of_interest %in% ff_channels)) {
        ff_markers <- markernames(ff)
        
        if (!all(markers_of_interest %in% ff_markers)) stop(paste0("Cannot find channels/markers of interest in ", file))
        
        channel_idcs       <- which(markers_of_interest %in% ff_markers)
        ff@exprs           <- ff@exprs[, channel_idcs]
        colnames(ff@exprs) <- channels_of_interest
      }
      
      classification      <- FlowSOM:::MapDataToCodes(som.grid$map$codes, ff@exprs)
      classification[, 1] <- som.metaClustering[classification[, 1]]
      
      freqs           <- table(classification[, 1])
      ordering        <- order(as.numeric(names(freqs)))
      freqs           <- freqs[ordering]
      freqs           <- rbind(as.numeric(names(freqs)), freqs)
      rownames(freqs) <- NULL; colnames(freqs) <- NULL
      absent_bins     <- which(!1:n_bins %in% freqs[1, ])
      if (length(absent_bins > 0)) {
        for (j in absent_bins) freqs <- cbind(freqs, c(j, 0))
        ordering <- order(freqs[1, ])
        freqs <- freqs[, ordering]
      }
      
      bins[[i]]$assignment <- classification[, 1]
      bins[[i]]$fp.abs     <- freqs[2, ]
      bins[[i]]$fp.rel     <- freqs[2, ] / sum(freqs[2, ])
    }
    names(bins) <- basename(fcs.paths)
    
    model                             <- som.grid
    model$metaclustering              <- as.factor(som.metaClustering)
    aggregate_matrix_binning          <- FlowSOM:::MapDataToCodes(som.grid$map$codes, aa_s3$agg$aggregate)
    aggregate_matrix_binning          <- som.metaClustering[aggregate_matrix_binning[, 1]]
    balanced_aggregate_matrix_binning <- som.metaClustering[som.grid$map$mapping[, 1]]
    
    rm(som.grid)
    
    
  } else if (method == "CellCnn") { ###############################
    
    ## Reload reticulate
    if ("reticulate" %in% .packages()) detach("package:reticulate", unload = TRUE, force = TRUE)    
    library(reticulate)
    
    cat("Method: CellCnn\n")
    
    outdir <- if (!is.null(method.settings$outdir)) { method.settings$outdir } else { "cellcnn_model" }
    
    ## Generate pre-processed .fcs files without junk cells (if needed)
    if (!is.null(junk)) {
      
      cleaned_fcs_folder <- if (!is.null(method.settings$cleaned_fcs_folder)) { method.settings$cleaned_fcs_folder } else { "cellcnn_fcs_files_cleaned" }
      
      if (!dir.exists(cleaned_fcs_folder)) dir.create(cleaned_fcs_folder)
      
      cat("Cleaning .fcs files...\n")
      for (path in fcs.paths) {
        file    <- basename(path)
        newpath <- file.path(cleaned_fcs_folder, file)
        ff      <- read.clean_FCS(path, junk = junk, id = file)
        
        cat(newpath, "\n")
        write.FCS(x = make_valid_fcs(ff@exprs, desc1 = ff@parameters$desc), filename = newpath)
      }
    }
    
    ## Create .csv input with names of input .fcs files and their labels
    fcs_csv_name      <- ".cellcnn_fcs_samples_with_labels.csv"
    
    fcs.paths.cellcnn <- if (!is.null(junk)) {
      sapply(fcs.paths, function(x) file.path(cleaned_fcs_folder, basename(x)))
    } else fcs.paths
    
    fcs_csv <- data.frame(fcs_filename = fcs.paths.cellcnn,
                          label = as.numeric(1:length(fcs.paths.cellcnn) %in% idcs.affected))
    write.csv(fcs_csv, fcs_csv_name, row.names = FALSE)
    cat("Input files and their labels:\n")
    print(fcs_csv)
    
    # Use correct Conda environment
    env_name <- if (!is.null(method.settings$env_name)) { method.settings$env_name } else { "tidycell_cellcnn_env" }
    use_condaenv(env_name)
    
    print(py_config())
    
    ## Run the analysis...
    
    nsamp_min <- min(length(idcs.controls), length(idcs.affected))
    
    cat("Preparing training data...\n")
    source_python(file.path(tidycell_path, "cellcnn_r_compatibles/prepare_data.py"))
    params.prepare_data <- list(fcs = fcs_csv_name,
                                marker_names = markers_of_interest,
                                quant_normed = FALSE,
                                arcsinh = FALSE,
                                cofactor = 5,
                                indir = "./",
                                seed = 1L,
                                train_perc = 0.75,
                                regression = FALSE,
                                per_sample = FALSE,
                                n_splits = if (nsamp_min < 4) { nsamp_min } else { NULL })
    d <- do.call(prepare_data, params.prepare_data)
    
    train_samples    <- d[[1]]
    train_phenotypes <- d[[2]]
    valid_samples    <- d[[3]]
    valid_phenotypes <- d[[4]]
    marker_names     <- d[[5]]
    fcs_info         <- d[[6]]
    samples          <- d[[7]]
    rm(d)
    
    for (i in 1:5) gc()
    
    scale             <- if (!is.null(method.settings$scale)) { method.settings$scale } else { TRUE }
    nfilter_choice    <- if (!is.null(method.settings$nfilter_choice)) { method.settings$nfilter_choice } else { 3L:9L }
    dendrogram_cutoff <- if (!is.null(method.settings$dendrogram_cutoff)) { method.settings$dendrogram_cutoff } else { 0.4 }
    max_epochs        <- if (!is.null(method.settings$max_epochs)) { method.settings$max_epochs } else { 20L }
    
    cat("NFILTER CHOICE:", nfilter_choice, "\n")
    
    cat("Creating model...\n")
    source_python(file.path(tidycell_path, "cellcnn_r_compatibles/model.py"))
    params.model <- list(ncell = 1000L,
                         nsubset = 1000L,
                         per_sample = TRUE,
                         subset_selection = "random",
                         scale = TRUE,
                         # quant_normed = FALSE,
                         # maxpool_percentages = c(0.01, 1, 5, 20, 100),
                         nfilter_choice = nfilter_choice,
                         # nrun = 15L,
                         # regression = FALSE,
                         # learning_rate = 0.005,
                         # coeff_l1 = 0,
                         # coeff_l2 = 0.0001,
                         max_epochs = max_epochs,
                         # patience = 5L,
                         dendrogram_cutoff = dendrogram_cutoff,
                         # accur_thres = 0.95,
                         verbose = 1L)
    model <- do.call(CellCnn, params.model)
    
    for (i in 1:5) gc()
    
    cat("Fitting to data...\n")
    model$fit(train_samples = train_samples,
              train_phenotypes = train_phenotypes,
              valid_samples = valid_samples,
              valid_phenotypes = valid_phenotypes,
              outdir = outdir)
    results <- model$results
    
    for (i in 1:5) gc()
    
    select_discriminative_filters <- if (!is.null(method.settings$select_discriminative_filters)) { method.settings$select_discriminative_filters } else { TRUE }
    
    cat("Extracting labels of aberrant populations...\n")
    source_python(file.path(tidycell_path, "cellcnn_r_compatibles/get_aberrant_pops.py"))
    params.get_aberrant_pops <- list(results_pointer = results,
                                     fcs = fcs_csv_name,
                                     samples_pointer = samples,
                                     marker_names = marker_names,
                                     filter_diff_thres = 0.2,
                                     filter_response_thres = 0,
                                     select_discriminative_filters = select_discriminative_filters,
                                     stat_test = NULL,
                                     group_a = "controls",
                                     group_b = "affected",
                                     group_names = NULL,
                                     regression = FALSE,
                                     train_samples = train_samples,
                                     train_phenotypes = train_phenotypes,
                                     valid_samples = valid_samples,
                                     valid_phenotypes = valid_phenotypes)
    aberrant_pops <- do.call(get_aberrant_pops, params.get_aberrant_pops)
    
    #save.image(file = "temp.RData")
    saveRDS(aberrant_pops, file = "aberrant_pops.RDS")
    
    cols <- colnames(aberrant_pops[[1]])
    cols.binary <- cols[grep("binary", cols)]
    cols.continuous <- cols[grep("continuous", cols)]
    
    responses <- do.call(rbind, aberrant_pops)[cols.continuous]
    
    if (!is.null(output_directory)) {
      ## Hierarchical clustering of filter responses
      colnames(responses) <- 1:ncol(responses)
      d <- distance_matrix(responses, dist_fun = "correlation")
      d <- 1 - d
      
      dendro <- as.dendrogram(hclust(d, "average"))
      
      png(file.path(output_directory, "filters_dendrogram.png"),
          width = 800,
          height = 800)
      layout(matrix(1))
      
      plot(dendro, main = "Clustering of CellCnn filters by pairwise correlations of filter responses", xlab = "Filter number", ylab = "1 - Pearson correlation")
      
      dev.off()
    }
    
    n_bins <- pop_count <- length(cols.binary)
    
    cat("Number of interesting filters:", n_bins, "\n")
    
    file.remove(fcs_csv_name)
    
    for (i in 1:length(aberrant_pops)) {
      binary     <- aberrant_pops[[i]][, cols.binary, drop = FALSE]
      continuous <- aberrant_pops[[i]][, cols.continuous, drop = FALSE]
      
      if (ncol(binary) > 1) {
        for (j in 2:ncol(binary)) binary[binary[, j] == 1, j] <- j
      }
      
      
      N     <- nrow(binary)
      freqs <- apply(binary, 2, function(x) sum(x > 0))
      
      bins[[i]]$fp.abs     <- freqs
      bins[[i]]$fp.rel     <- freqs / N
      bins[[i]]$assignment <- binary
      bins[[i]]$continuous <- continuous
    }
    names(bins) <- basename(fcs.paths)
    
    bin_overlap <- TRUE # 'bins' can possibly overlap
    
    aggregate_matrix_binning <- lapply(1:length(bins), function(idx) {
      colnum <- n_bins
      
      if (is.null(junk)) {
        assignment <- bins[[idx]]$assignment
      } else {
        assignment                      <- matrix(0, nrow = nrow(read.FCS(fcs.paths[idx])), ncol = colnum)
        j                               <- unique(unlist(junk[[basename(fcs.paths[idx])]]))
        which.not_removed               <- 1:nrow(assignment)
        if (length(j) > 0) which.not_removed <- which.not_removed[-j]
        assignment[which.not_removed, ] <- as.matrix(bins[[idx]]$assignment)
      }
      assignment[aa_s3$agg$downsample_idcs.aggregate[[idx]], ]
    }) %>% do.call(rbind, .)
    
    controls_matrix_binning <- lapply(1:length(idcs.controls), function(i) {
      idx <- idcs.controls[i]
      colnum <- n_bins
      if (is.null(junk)) {
        assignment <- bins[[idx]]$assignment
      } else {
        assignment                      <- matrix(0, nrow = nrow(read.FCS(fcs.paths[idx])), ncol = colnum)
        j                               <- unique(unlist(junk[[basename(fcs.paths[idx])]]))
        which.not_removed               <- 1:nrow(assignment)
        if (length(j) > 0) which.not_removed <- which.not_removed[-j]
        assignment[which.not_removed, ] <- as.matrix(bins[[idx]]$assignment)
      }
      assignment[aa_s3$agg$downsample_idcs.controls[[i]], ]
    }) %>% do.call(rbind, .)
    
    affected_matrix_binning <- lapply(1:length(idcs.affected), function(i) {
      idx <- idcs.affected[i]
      colnum <- n_bins
      if (is.null(junk)) {
        assignment <- bins[[idx]]$assignment
      } else {
        assignment                      <- matrix(0, nrow = nrow(read.FCS(fcs.paths[idx])), ncol = colnum)
        j                               <- unique(unlist(junk[[basename(fcs.paths[idx])]]))
        which.not_removed               <- 1:nrow(assignment)
        if (length(j) > 0) which.not_removed <- which.not_removed[-j]
        assignment[which.not_removed, ] <- as.matrix(bins[[idx]]$assignment)
      }
      assignment[aa_s3$agg$downsample_idcs.affected[[i]], ]
    }) %>% do.call(rbind, .)
    
    balanced_aggregate_matrix_binning <- rbind(controls_matrix_binning, affected_matrix_binning)
    
    model <- list(n_bins = n_bins)
    
  } else {
    stop("Invalid method")
  }
  
  ## Compute aberrances based on controls fingerprint baseline
  baseline <- lapply(bins[idcs.controls], function(x) x$fp.rel) %>% do.call(rbind, .)
  baseline <- apply(baseline, 2, median)
  
  aberrance <- vector(mode = "list", length = length(fcs.paths)); names(aberrance) <- basename(fcs.paths)
  for (i in 1:length(fcs.paths)) {
    fp                  <- bins[[i]]$fp.rel
    aberrance[[i]]$real <- fp - baseline
    aberrance[[i]]$plot <- .sgn(fp - baseline) * sqrt(abs(fp - baseline))
  }
  names(aberrance) <- basename(fcs.paths)
  
  vals.controls <- lapply(aberrance[idcs.controls], function(x) x$real)
  vals.affected <- lapply(aberrance[idcs.affected], function(x) x$real)
  
  d0 <- lapply(vals.controls, function(x) data.frame(bin = as.factor(1:length(x)), aberrance = x))
  d1 <- lapply(vals.affected, function(x) data.frame(bin = as.factor(1:length(x)), aberrance = x))
  d0 <- do.call(rbind, d0)
  d1 <- do.call(rbind, d1)
  
  p.vals <- vector(mode = "numeric", length = length(unique(d0$bin)))
  for (i in 1:length(p.vals)) {
    vals0 <- d0$aberrance[d0$bin == i]
    vals1 <- d1$aberrance[d1$bin == i]
    suppressWarnings(p.vals[i] <- wilcox.test(vals0, vals1, alternative = "two.sided")$p.value)
  }
  
  if (!bin_overlap) {
    MFIs <- lapply(sort(unique(balanced_aggregate_matrix_binning)), function(bin) {
      bin_idcs <- which(balanced_aggregate_matrix_binning == bin)
      if (length(bin_idcs) > 1) {

        submatrix <- aa_s3$agg$balanced_aggregate[bin_idcs, channels_of_interest]

        return(apply(submatrix, 2, median))
      } else {
        return(NULL)
      }
    }) %>% do.call(rbind, .)

    colnames(MFIs) <- markers_of_interest
    rownames(MFIs) <- 1:nrow(MFIs)
    MFIs <- t(MFIs)
    
  } else { ## IF bin_overlap == TRUE:
    MFIs <- lapply(1:ncol(balanced_aggregate_matrix_binning), function(bin) {
      bin_idcs <- which(balanced_aggregate_matrix_binning[, bin] > 0)
      if (length(bin_idcs) > 1) {
        
        submatrix <- aa_s3$agg$balanced_aggregate[bin_idcs, channels_of_interest]
        
        return(apply(submatrix, 2, median))
      } else {
        return(NULL)
      }
    }) %>% do.call(rbind, .)
    
    colnames(MFIs) <- markers_of_interest
    rownames(MFIs) <- 1:nrow(MFIs)
    MFIs <- t(MFIs)
  }
  
  list(method = method,
       bin_overlap = bin_overlap,
       aggregate_matrix_binning = aggregate_matrix_binning,
       balanced_aggregate_matrix_binning = balanced_aggregate_matrix_binning,
       MFIs = MFIs,
       p.vals = p.vals,
       fcs.paths.original = fcs.paths.original,
       fcs.paths = fcs.paths,
       idcs.controls.original = idcs.controls.original,
       idcs.controls = idcs.controls,
       idcs.affected.original = idcs.affected.original,
       idcs.affected = idcs.affected,
       idcs.unused = idcs.unused,
       model = model,
       n_bins = n_bins,
       binning = bins,
       aberrance = aberrance)
}

## The MEM (Marker Enrichment Modelling) equation comes from:
### Diggins. K. E., Gandelman, J. S., Roe, C. E., & Irish, J. M. (2018). Generating quantitative cell identity labels
### with marker enrichment modeling (MEM). Current Protocols in Cytometry, 83, 10.21.1-10.21.28. doi: 10.1002/cpcy.34
.mem <- function(MAGpop, MAGref, IQRpop, IQRref) (abs(MAGpop - MAGref) + IQRref/IQRpop - 1) * .sgn(MAGpop - MAGref)

## Function: provided .fcs data from controls & affected, compute MEM scores for each marker of interest.
## The MEM (Marker Enrichment Modelling) equation comes from:
### Diggins. K. E., Gandelman, J. S., Roe, C. E., & Irish, J. M. (2018). Generating quantitative cell identity labels
### with marker enrichment modeling (MEM). Current Protocols in Cytometry, 83, 10.21.1-10.21.28. doi: 10.1002/cpcy.34
mem.scores <- function(aa_s3,
                       fcs.paths,
                       idcs.unused = NULL,
                       junk = NULL,
                       idx = NULL,
                       channels_of_interest,
                       markers_of_interest) {
  
  cat("Computing MEM scores for each marker\n")
  
  efcs.controls            <- aa_s3$agg$controls_matrix
  efcs.controls            <- efcs.controls[, channels_of_interest]
  colnames(efcs.controls)  <- markers_of_interest
  
  efcs.affected            <- aa_s3$agg$affected_matrix
  efcs.affected            <- efcs.affected[, channels_of_interest]
  colnames(efcs.affected)  <- markers_of_interest
  
  efcs.aggregate           <- aa_s3$agg$balanced_aggregate
  efcs.aggregate           <- efcs.aggregate[, channels_of_interest]
  colnames(efcs.aggregate) <- markers_of_interest
  
  median.controls          <- apply(efcs.controls, 2, median)
  median.affected          <- apply(efcs.affected, 2, median)
  median.aggregate         <- apply(efcs.aggregate, 2, median)
  iqr.controls             <- apply(efcs.controls, 2, IQR)
  iqr.affected             <- apply(efcs.controls, 2, IQR)
  iqr.aggregate            <- apply(efcs.aggregate, 2, IQR)
  
  scores.aggregated <- .mem(median.affected, median.controls, iqr.affected, iqr.controls)
  
  if (!is.null(idcs.unused)) fcs.paths <- fcs.paths[-idcs.unused]
  
  scores.per_sample <- vector(mode = "list", length = length(fcs.paths)); names(scores.per_sample) <- fcs.paths
  for (i in 1:length(fcs.paths)) {
    file           <- fcs.paths[i]
    ff             <- read.clean_FCS(file, junk = junk, id = basename(file))
    efcs           <- ff@exprs
    efcs           <- efcs[, channels_of_interest]
    colnames(efcs) <- markers_of_interest
    
    scores.this_sample <- vector(mode = "list", length = length(markers_of_interest)); names(scores.this_sample) <- markers_of_interest
    for (marker in markers_of_interest) {
      marker_values      <- efcs[, marker]
      median.this_sample <- median(marker_values)
      iqr.this_sample    <- IQR(marker_values)
      
      scores.this_sample[[marker]] <- .mem(median.this_sample, median.aggregate[marker], iqr.this_sample, iqr.aggregate[marker])
    }
    scores.this_sample        <- unlist(scores.this_sample)
    names(scores.this_sample) <- markers_of_interest
    scores.per_sample[[i]]    <- scores.this_sample
  }
  
  list(scores.aggregated = scores.aggregated,
       scores.per_sample = scores.per_sample)
}
