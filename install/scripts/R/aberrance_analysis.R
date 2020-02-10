## aberrance analysis module

## Function: run aberrance analysis on specified .fcs files with given parameters
aberrance_analysis <- function(fcs.directory,
                               fcs.files,
                               indices.controls,
                               indices.affected,
                               indices.unused      = NULL,
                               markers.analysis,
                               channels.analysis,
                               additional_columns  = NULL,              # annotation columns
                               flowjo_workspace    = NULL,              # path
                               exclusive_gates     = NULL,              # desired substrings for filtering nodes in gating hierarchy
                               inv_gate_indices    = NULL,
                               binning_methods     = list(list(method   = "hyperbinning",
                                                               n_bins   = 32),
                                                          list(method   = "FlowSOM",
                                                               n_bins   = 32,
                                                               settings = list(grid.x = 15,
                                                                               grid.y = 15,
                                                                               scale  = FALSE)),
                                                          list(method   = "hyperbinning",
                                                               n_bins   = 512),
                                                          list(method   = "CellCnn",
                                                               settings = list(env_name       = "tidycell_cellcnn_env",
                                                                               outdir         = "cellcnn_model",
                                                                               nfilter_choice = 14L:15L))),
                               tidycell_path       = "./scipts/",       # path to R and Python scripts
                               output_directory,
                               concatenated_fcs    = NULL,              # name
                               precomputed_dimred  = NULL,
                               precomputed_concat  = NULL,
                               dimension_reduction = list(list(method   = "umap",
                                                               markers  = markers.analysis,
                                                               channels = channels.analysis)),
                               visual_sanity_checks         = FALSE,
                               files_preprocessed_already   = FALSE,
                               asinh_transform              = FALSE,
                               asinh_cofactor               = 1/5,
                               estimateLogicle_transform = estimateLogicle_transform,
                               preprocessed_files_directory = ".fcs_files_preprocessed",
                               seed = 1) {
  
  set.seed(seed)
  
  if (!is.null(inv_gate_indices)) {
    o <- sapply(fcs.files, function(x) which(names(inv_gate_indices) == x))
    if (is.list(o) || length(o) != length(fcs.files)) {
      if (length(o) == length(fcs.files)) {
        message('inv_gate_indices incongruent, using indices')
        names(inv_gate_indices) <- fcs.files
      } else {
        stop('inv_gate_indices length incongruent')
      }
    } else {
      inv_gate_indices <- inv_gate_indices[o]
    }
  }
  
  cat("Controls samples indices:", paste(indices.controls, collapse = ", "), "\n")
  cat("Affected samples indices:", paste(indices.affected, collapse = ", "), "\n")
  cat("Unused samples indices:",   paste(indices.unused,   collapse = ", "), "\n")
  
  if (length(indices.unused) == 1 && indices.unused == 0) indices.unused <- NULL
  
  if (!file.exists(preprocessed_files_directory)) dir.create(preprocessed_files_directory)
  if (!file.exists(output_directory))             dir.create(output_directory)
  
  fcs.paths <- file.path(fcs.directory, fcs.files)
  
  ## Mark pre-gated events (if applicable)
  if (is.null(inv_gate_indices)) {
    not_gated <- NULL
    if (!is.null(exclusive_gates)) {
      not_gated <- get_gate_indices(fcs.dir    = fcs.directory,
                                    fcs.files  = fcs.files,
                                    wsp        = flowjo_workspace,
                                    gate_names = exclusive_gates,
                                    invert     = TRUE)
    }
  } else {
    not_gated <- inv_gate_indices
  }
  
  ## Transform data & remove irrelevant events (based on pre-gated populations,
  ## preprocess_fcs can be called separately for QC, removing extrema, etc.).
  junk <- 0
  fcs.paths.analysis <- fcs.paths
  if (asinh_transform || estimateLogicle_transform || !is.null(to_remove) || !is.null(not_gated)) {
    if (!files_preprocessed_already) {
      channels     <- colnames(read.FCS(fcs.paths.analysis[1]))
      idcs.scatter <- c(grep("FSC-", channels), grep("SSC-", channels))
      if (length(idcs.scatter) > 0) {
        message(paste0('Identified scatter channels ', paste(channels[idcs.scatter], collapse = ', ')));
        channels.transform <- channels[-idcs.scatter]  
      } else {
        channels.transform <- channels
      }
      channels.transform <- channels.transform[!channels.transform %in% c("TIME", "Time")]
      transformation <- if (asinh_transform) {
        transformList(channels.transform, arcsinhTransform(a = 0, b = asinh_cofactor, c = 0))
      } else if (estimateLogicle_transform) {
        'estimateLogicle'
      } else { NULL }
      preprocess_fcs(fcs.paths,
                     scatter_idcs         = if (length(idcs.scatter) > 0) { idcs.scatter } else { NULL },
                     known_junk           = not_gated,
                     transformation       = transformation,
                     foldername           = preprocessed_files_directory,
                     out.events_to_remove = junk)
    } else { # IF files_preprocessed_already == TRUE:
      junk <- readRDS('.junk.RDS')
    }
    fcs.paths.analysis <- file.path(preprocessed_files_directory, fcs.files)
    
    saveRDS(junk, file = file.path(output_directory, paste0(output_directory, '_junk_events_indices.RDS')))
  }
  
  
  if (length(junk) == 0 || sum(sapply(junk, length)) == 0) junk <- NULL
  saveRDS(junk, ".junk.RDS")
  
  ## Create an environment for storing large data (and avoid passing them by copy)
  aa_s3 <- aa()
  
  ## Concatenate data
  if (is.null(precomputed_concat)) {
    aa_s3$concat <- concatenate_fcs(fcs.paths.analysis, junk = junk)
    saveRDS(aa_s3$concat, ".concat.RDS")
  } else {
    aa_s3$concat <- precomputed_concat
    rm(precomputed_concat)
  }
  
  ## Dimension reduction
  if (is.null(precomputed_dimred)) {
    aa_s3$dimred.list <- lapply(dimension_reduction, function(d) {
      dimred.layout(aa_s3                = aa_s3,
                    channels_of_interest = d$channels,
                    method               = d$method,
                    method.settings      = d$method.settings,
                    normalise.quantile   = .9)
    })
  } else {
    aa_s3$dimred.list <- precomputed_dimred
    rm(precomputed_dimred)
  }
  if (!is.null(aa_s3$dimred.list)) {
    saveRDS(aa_s3$dimred.list, file = file.path(output_directory, paste0(output_directory, "_dimred.RDS")))
    saveRDS(aa_s3$dimred.list, file = ".dimred.RDS")
  }
  
  ## Stratified downsampling for model building
  aa_s3$agg <- aggregate_data(fcs.paths.analysis,
                              junk          = junk,
                              idcs.controls = indices.controls,
                              idcs.affected = indices.affected,
                              idcs.unused   = indices.unused)
  
  ## Identify aberrant populations and compute their relative over- or underrepresentation
  aberrances <- lapply(binning_methods, function(bm) {
    aberrance(aa_s3                = aa_s3,
              junk                 = junk,
              fcs.paths            = fcs.paths.analysis,
              idcs.controls        = indices.controls,
              idcs.affected        = indices.affected,
              idcs.unused          = indices.unused,
              channels_of_interest = channels.analysis,
              markers_of_interest  = markers.analysis,
              method               = bm$method,
              method.settings      = bm$settings,
              tidycell_path        = tidycell_path,
              n_bins               = bm$n_bins,
              output_directory     = output_directory)
  })
  
  ## Save aberrances
  saveRDS(aberrances, file = file.path(output_directory, paste0(output_directory, "_aberrances.RDS")))
  
  names(aberrances) <- lapply(1:length(binning_methods), function(idx) {
    bm <- binning_methods[[idx]]
    if (!is.null(bm$n_bins)) {
      return(paste0(idx, "_", bm$method, bm$n_bins))
    } else {
      return(paste0(idx, "_", bm$method))
    }
  })
  
  ## Compute MEM scores (measures of marker enrichment)
  mem_scores <- mem.scores(aa_s3                = aa_s3,
                           fcs.paths            = fcs.paths.analysis,
                           junk                 = junk,
                           idcs.unused          = indices.unused,
                           channels_of_interest = channels.analysis,
                           markers_of_interest  = markers.analysis)
  
  ## Write enhanced .fcs files with additional columns (main output of this analysis)
  fcs.write_enhanced(fcs.paths                  = fcs.paths,
                     aberrances.list            = aberrances,
                     concatenated_fcs           = concatenated_fcs,
                     additional_columns         = additional_columns,
                     dimred.list                = aa_s3$dimred.list,
                     junk                       = junk,
                     removed_flag               = -1000,
                     clean_indices_for_checking = aa_s3$concat$sample_idcs,
                     output.path                = output_directory)
  
  cat("Done enhancing .fcs files\n")
  
  ## Create visual outputs (if enabled)
  if (visual_sanity_checks) {
    cat("Generating visuals\n")
    visuals_path <- file.path(output_directory, "sanity_checks")
    if (!file.exists(visuals_path)) { dir.create(visuals_path) }
    
    cat("Plotting MEM scores \n")
    create_visual.mem_scores(scores        = mem_scores,
                             idcs.controls = indices.controls,
                             idcs.affected = indices.affected,
                             idcs.unused   = indices.unused,
                             visual.path   = file.path(visuals_path, "mem_scores.png"))
    
    cat("Plotting aberrance box plots, bin size bar plots and bin expression heatmaps\n")
    for (idx in 1:length(aberrances)) {
      ab <- aberrances[[idx]]
      if (ab$method != "CellCnn") {
        ab   <- aberrances[[idx]]
        name <- paste0(idx, "_", ab$method, "_", ab$n_bins)
        bin_order <- 0
        create_visual.bin_heatmap(aberrances    = ab,
                                  out.bin_order = bin_order,
                                  visual.path   = file.path(visuals_path, paste0(name, "_heatmap.png")))
        create_visual.bin_sizes(aberrances  = ab,
                                title       = paste0(name, " bin sizes"),
                                bin.order   = bin_order,
                                visual.path = file.path(visuals_path, paste0(name, "_binsizes.png")))
        create_visual.aberrance(aberrances = ab,
                                title      = paste0(name, " aberrances (alpha = 0.1)"),
                                #bin.order  = bin_order,
                                alpha      = 0.1,
                                visual.path   = file.path(visuals_path, paste0(name, "_aberrances.png")))
      }
    }
  }
}

## Function: parse .csv inputs and and pass them as parameters to aberrance analysis function
call.aberrance_analysis <- function(csv.analysis_inputs,
                                    csv.analysis_markers,
                                    csv.analysis_settings,
                                    tidycell_path        = "./",
                                    inv_gate_indices     = NULL) {
  
  if (!file.exists(csv.analysis_inputs))   stop(paste0('File ', csv.analysis_inputs, ' not found. Terminating analysis'))
  if (!file.exists(csv.analysis_markers))  stop(paste0('File ', csv.analysis_markers, ' not found. Terminating analysis'))
  if (!file.exists(csv.analysis_settings)) stop(paste0('File ', csv.analysis_settings, ' not found. Terminating analysis'))
  
  ## Parse analysis inputs (classify .fcs files for aberrance analysis/dimension reductions,
  ## parse information to add to enhanced .fcs files as additional columns (with constant value per file)
  inputs           <- read.table(csv.analysis_inputs, sep = ",", header = TRUE)
  which_cols.class <- grep("^class", colnames(inputs))
  
  if (length(which_cols.class) < 1) stop('Analysis inputs table must contain "class" columns. Terminating analysis')
  
  indices.controls <- vector(mode = "list", length = length(which_cols.class))
  indices.affected <- vector(mode = "list", length = length(which_cols.class))
  indices.unused   <- vector(mode = "list", length = length(which_cols.class))
  
  fcs.files <- vector(mode = "list", length = length(which_cols.class))
  for (idx.class in 1:length(which_cols.class)) {
    na                            <- which(is.na(inputs[, which_cols.class[idx.class]]) | inputs[, which_cols.class[idx.class]] == 'NA')
    indices.unused[[idx.class]]   <- if (length(na) == 0) { 0 } else { na }
    
    inp                           <- inputs[, which_cols.class[idx.class]]
    indices.controls[[idx.class]] <- which(inp == 0)
    indices.affected[[idx.class]] <- which(inp == 1)
    
    if (any(!na.omit(inp) %in% c(0, 1))) stop('Class labels in analysis inputs table must be 0 and 1. Terminating analysis')
    
    fcs.files[[idx.class]] <- as.character(inputs$fcs_file)
  }
  
  additional_columns <- NULL
  for (idx in which_cols.class) {
    name                       <- colnames(inputs)[idx]
    vals                       <- inputs[[name]]
    vals[is.na(vals)]          <- -1
    additional_columns[[name]] <- vals
  }
  if (max(which_cols.class) < ncol(inputs)) {
    which_cols.additional <- max(which_cols.class + 1):ncol(inputs)
    for (idx in which_cols.additional) additional_columns[[colnames(inputs)[idx]]] <- as.character(inputs[, idx])
  }
  
  ## Parse markers (& channels for analysis)
  markers           <- read.table(csv.analysis_markers, sep = ",", header = TRUE)
  markers.analysis  <- as.character(markers$markers[markers$analysis == 1])
  if (length(markers.analysis) < 1) stop('No markers are selected for analysis in the markers table. Terminating analysis')
  channels.analysis <- as.character(markers$channels[markers$analysis == 1])
  markers.tsne      <- as.character(markers$markers[markers$tsne == 1])
  channels.tsne     <- as.character(markers$channels[markers$tsne == 1])
  markers.umap      <- as.character(markers$markers[markers$umap == 1])
  channels.umap     <- as.character(markers$channels[markers$umap == 1])
  markers.embedsom  <- as.character(markers$markers[markers$embedsom == 1])
  channels.embedsom <- as.character(markers$channels[markers$embedsom == 1])
  
  settings <- read.table(csv.analysis_settings, sep = ",", header = TRUE, na.strings = "")
  
  cn <- paste(colnames(settings), collapse = ', ')
  cn_correct <- 'fcs_directory, workspace, exclusive_gates, transform, transform_cofactor, binning_methods, output_directory, concatenate_fcs, dimension_reduction, precomputed_dimension_reduction, sanity_checks'
  if (cn != cn_correct)
    stop('Invalid analysis_settings columns. Please check format of your analysis settings table. Terminating analysis')
  
  .clean  <- function(vec) as.character(na.omit(vec))
  .method <- function(str) {
    n_bins <- 32
    input  <- strsplit(str, "=")[[1]]
    if (length(input) == 2) n_bins <- input[2]
    
    n_bins <- as.numeric(as.character(n_bins))
    
    list(method = toupper(input[1]),
         n_bins = as.numeric(n_bins))
  }
  
  fcs.directory    <- .clean(settings$fcs_directory)
  flowjo_workspace <- .clean(settings$workspace)
  
  if (length(flowjo_workspace) == 0) flowjo_workspace <- NULL
  
  exclusive_gates <- .clean(settings$exclusive_gates)
  if (length(exclusive_gates) == 0) exclusive_gates <- NULL
  
  asinh_transform           <- FALSE
  estimateLogicle_transform <- FALSE
  tt <- .clean(settings$transform)
  if (length(tt) > 0) {
    asinh_transform           <- (tt == "asinh" || tt == "arcsinh")
    estimateLogicle_transform <- (tt == "estimateLogicle")
  } else {
    stop('Invalid transformation entered in the analysis settings table. Terminating analysis')
  }
  asinh_cofactor <- as.numeric(.clean(settings$transform_cofactor))
  
  methods               <- lapply(.clean(settings$binning_methods), .method)
  hyperbinning.defaults <- NULL
  FlowSOM.defaults      <- list(grid.x = 15,
                                grid.y = 15,
                                scale = FALSE)
  CellCnn.defaults      <- list(outdir = "cellcnn_model",
                                nfilter_choice = 3L:15L,
                                select_discriminative_filters = FALSE,
                                max_epochs = 20L)
  
  binning_methods <- lapply(methods, function(method) {
    if (method$method == "HYPERBINNING") {
      return(list(method   = "hyperbinning",
                  settings = hyperbinning.defaults,
                  n_bins   = method$n_bins))
    } else if (method$method == "FLOWSOM") {
      return(list(method   = "FlowSOM",
                  settings = FlowSOM.defaults,
                  n_bins   = method$n_bins))
    } else if (method$method == "CELLCNN") {
      return(list(method   = "CellCnn",
                  settings   = CellCnn.defaults))
    } else {
      stop(paste0("Invalid binning method '", method$method, "' entered in the analysis settings table. Terminating analysis"))
    }
  })
  
  precomp_dimred <- .clean(settings$precomputed_dimension_reduction)
  precomp_dimred <- if (length(precomp_dimred) == 0) { '' } else { precomp_dimred[1] }
  
  output_directory <- .clean(settings$output_directory)
  concatenated_fcs <- if (.clean(settings$concatenate_fcs) == 1) { 'concat.fcs' } else { NULL }
  
  dimreds <- .clean(settings$dimension_reduction)
  dimension_reduction <- lapply(dimreds, function(dimred) {
    if (toupper(dimred) == "TSNE") {
      return(list(method   = "tsne",
                  markers  = markers.tsne,
                  channels = channels.tsne))
    } else if (toupper(dimred) == "UMAP") {
      return(list(method   = "umap",
                  markers  = markers.umap,
                  channels = channels.umap,
                  method.settings = list(n_neighbours = 50L)))
    } else if (toupper(dimred) == "EMBEDSOM" || toupper(dimred) == "EMBEDSOM:FLAT") {
      return(list(method   = "embedsom",
                  markers  = markers.embedsom,
                  channels = channels.embedsom,
                  method.settings = list(grid.x   = 20,
                                         grid.y   = 20,
                                         emcoords = "flat")))
    } else if (toupper(dimred) == "EMBEDSOM:UMAP") {
      return(list(method   = "embedsom",
                  markers  = markers.embedsom,
                  channels = channels.embedsom,
                  method.settings = list(grid.x  = 20,
                                         grid.y  = 20,
                                         emcoords = "uwot::umap")))
    } else if (toupper(dimred) == "EMBEDSOM:TSNE") {
      return(list(method   = "embedsom",
                  markers  = markers.embedsom,
                  channels = channels.embedsom,
                  method.settings = list(grid.x   = 20,
                                         grid.y   = 20,
                                         emcoords = "tsne")))
    } else {
      stop(paste0("Invalid dimension reduction method '", dimred, "' entered in the analysis settings table. Terminating analysis"))
    }
  })
  
  fileConn <- file("input_check_result.txt")
  error    <- input_check(csv.analysis_inputs  = csv.analysis_inputs,
                          csv.analysis_markers = csv.analysis_markers,
                          fcs.directory        = fcs.directory,
                          file_connection      = fileConn)
  close(fileConn)
  if (error) {                  
    stop("Please check the pairing of channels and markers in the markers table and your .fcs files. For more info on mismatches, read through 'input_check_result.txt' Terminating analysis")
  }
  
  visual_sanity_checks <- .clean(settings$sanity_checks) == '1'
  
  ## Compute required dimension reduction(s) in first iteration of analysis
  
  if (precomp_dimred != '') {
    if (!file.exists(precomp_dimred)) stop(paste0("Precomputed dimension reduction file '", precomp_dimred, "' not found. Terminating analysis"))
    message('Loading precomputed dimension reduction layout...')
    precomp_dimred <- readRDS(precomp_dimred)
    if (is.list(precomp_dimred) && length(precomp_dimred) == 0) stop('Precomputed dimension reduction object is an empty list. Terminating analysis')
  } else {
    precomp_dimred <- NULL
  }
  
  out <- paste0(output_directory, "_1")
  aberrance_analysis(fcs.directory              = fcs.directory,
                     fcs.files                  = fcs.files[[1]],
                     indices.controls           = indices.controls[[1]],
                     indices.affected           = indices.affected[[1]],
                     indices.unused             = indices.unused[[1]],
                     markers.analysis           = markers.analysis,
                     channels.analysis          = channels.analysis,
                     additional_columns         = additional_columns,
                     flowjo_workspace           = flowjo_workspace,
                     exclusive_gates            = exclusive_gates,
                     inv_gate_indices           = if (!is.null(inv_gate_indices)) { readRDS(inv_gate_indices) } else { NULL },
                     binning_methods            = binning_methods,
                     tidycell_path              = tidycell_path,
                     output_directory           = out,
                     concatenated_fcs           = paste0("1_", concatenated_fcs),
                     precomputed_dimred         = precomp_dimred,
                     precomputed_concat         = NULL,
                     dimension_reduction        = dimension_reduction,
                     visual_sanity_checks       = visual_sanity_checks,
                     files_preprocessed_already = FALSE,
                     asinh_transform            = asinh_transform,
                     asinh_cofactor             = asinh_cofactor,
                     estimateLogicle_transform  = estimateLogicle_transform)
  
  write.csv(inputs,   file.path(out, "analysis_inputs.csv"),   row.names = FALSE)
  write.csv(markers,  file.path(out, "analysis_markers.csv"),  row.names = FALSE)
  write.csv(settings, file.path(out, "analysis_settings.csv"), row.names = FALSE)
  
  if (length(which_cols.class) > 1) {
    if (is.null(precomp_dimred)) {
      if (file.exists(".dimred.RDS")) {
        precomp_dimred <- readRDS(".dimred.RDS")
      }
    }
    
    fileConn <- file(file.path(out, "analysis_id.txt"))
    writeLines(c("ABERRANCE ANALYSIS RESULTS", paste0("This is analysis corresponding to the classification column 1 in analysis_inputs.")), fileConn)
    close(fileConn)
    
    for (idx.class in 2:length(which_cols.class)) {
      out <- paste0(output_directory, "_", idx.class)
      aberrance_analysis(fcs.directory        = fcs.directory,
                         fcs.files            = fcs.files[[idx.class]],
                         indices.controls     = indices.controls[[idx.class]],
                         indices.affected     = indices.affected[[idx.class]],
                         indices.unused       = indices.unused[[idx.class]],
                         markers.analysis     = markers.analysis,
                         channels.analysis    = channels.analysis,
                         additional_columns   = additional_columns,
                         flowjo_workspace     = flowjo_workspace,
                         exclusive_gates      = exclusive_gates,
                         inv_gate_indices     = if (!is.null(inv_gate_indices)) { readRDS(inv_gate_indices) } else { NULL },
                         binning_methods      = binning_methods,
                         tidycell_path        = tidycell_path,
                         output_directory     = out,
                         concatenated_fcs     = paste0(idx.class, "_", concatenated_fcs),
                         precomputed_dimred   = precomp_dimred,
                         precomputed_concat   = readRDS(".concat.RDS"),
                         dimension_reduction  = dimension_reduction,
                         visual_sanity_checks = visual_sanity_checks,
                         files_preprocessed_already = TRUE,
                         asinh_transform      = asinh_transform,
                         asinh_cofactor       = asinh_cofactor,
                         estimateLogicle_transform = estimateLogicle_transform)
      write.csv(inputs,   file.path(out, "analysis_inputs.csv"),   row.names = FALSE)
      write.csv(markers,  file.path(out, "analysis_markers.csv"),  row.names = FALSE)
      write.csv(settings, file.path(out, "analysis_settings.csv"), row.names = FALSE)
      
      fileConn <- file(file.path(out, "analysis_id.txt"))
      writeLines(c("ABERRANCE ANALYSIS RESULTS", paste0("This is analysis corresponding to the classification column ", idx.class, " in analysis_inputs.")), fileConn)
      close(fileConn)
    }
  }
}
