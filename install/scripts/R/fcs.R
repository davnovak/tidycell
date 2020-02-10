## fcs enhancement module

get_fcs_paths <- function(dir) list.files(dir, full.names = TRUE, pattern = ".*fcs$", recursive = TRUE)

fcs.add_col <- function(ff, new_col, colname = "label") {
  
  efcs <- ff@exprs
  
  N <- nrow(efcs)
  len <- length(new_col)
  
  if (N != len)  stop(paste0("Number of rows of expression matrix is ", N, ", whereas length of new column is ", len, "."))
  
  params <- ff@parameters
  pd     <- pData(params)
  cols   <- as.vector(pd$name)
  idcs   <- match(cols, pd$name)
  
  if (any(is.na(idcs))) stop("Invalid column specifier")
  
  channel_number     <- ncol(ff) + 1
  channel_id         <- paste0("$P", channel_number)
  channel_name       <- colname
  channel_range      <- max(new_col) + 1
  channel_min        <- min(0, min(new_col) - 1)
  plist              <- matrix(c(channel_name, channel_name, channel_range,
                                 channel_min, channel_range - 1))
  rownames(plist)    <- c("name", "desc", "range", "minRange", "maxRange")
  colnames(plist)    <- c(channel_id)
  pd                 <- rbind(pd, t(plist))
  pData(params)      <- pd
  channel_names      <- colnames(efcs)
  efcs.mod           <- cbind(efcs, new_col)
  colnames(efcs.mod) <- c(channel_names, colname)
  
  ff.mod             <- flowFrame(efcs.mod, params, description = description(ff))
  
  keyval                                      <- list()
  keyval[[paste0("$P", channel_number, "B")]] <- "32"
  keyval[[paste0("$P", channel_number, "R")]] <- toString(channel_range)
  keyval[[paste0("$P", channel_number, "E")]] <- "0,0"
  keyval[[paste0("$P", channel_number, "N")]] <- channel_name
  keyval[[paste0("$P", channel_number, "S")]] <- channel_name
  keyword(ff.mod)                             <- keyval
  
  flowCoreP_Rmax <- paste0("flowCore_$P", channel_number, "Rmax")
  flowCoreP_Rmin <- paste0("flowCore_$P", channel_number, "Rmin")
  
  description(ff.mod)[flowCoreP_Rmax] <- max(20000, description(ff.mod)$`flowCore_$P1Rmax`)
  description(ff.mod)[flowCoreP_Rmin] <- 0
  
  return(ff.mod)
}

fcs.write_enhanced <- function(fcs.paths,
                               aberrances.list,
                               concatenated_fcs           = NULL,
                               additional_columns         = NULL,
                               dimred.list                = NULL,
                               junk                       = NULL,
                               removed_flag               =  -10,
                               clean_indices_for_checking = NULL,
                               output.path                = "fcs_files_enhanced") {
  
  ## Get sample indices for concatenated expression matrix
  sample_idcs       <- list()
  sample_idcs$full  <- vector(mode = "list", length = length(fcs.paths))
  sample_idcs$clean <- vector(mode = "list", length = length(fcs.paths))
  for (idx in 1:length(fcs.paths)) {
    path <- fcs.paths[idx]
    ff   <- read.FCS(path)
    if (idx == 1) {
      sample_idcs$full[[idx]]     <- c(1, nrow(ff@exprs))
      sample_idcs$clean[[idx]]    <- c(1, nrow(ff@exprs) - length(unlist(junk[[idx]])))
    } else {
      previous_end                <- sample_idcs$full[[idx - 1]][2]
      sample_idcs$full[[idx]]     <- c(previous_end + 1, previous_end + nrow(ff@exprs))
      
      previous_end <- sample_idcs$clean[[idx - 1]][2]
      sample_idcs$clean[[idx]]    <- c(previous_end + 1, previous_end + nrow(ff@exprs) - length(unlist(junk[[idx]])))
    }
    names(sample_idcs$clean)[idx] <- names(sample_idcs$full)[idx] <- basename(path)
  }
  
  if (!is.null(clean_indices_for_checking)) {
    for (i in 1:length(clean_indices_for_checking)) {
      if (!identical(clean_indices_for_checking[[i]], sample_idcs$clean[[i]])) stop("Clean event indices misaligned")
    }
  }
  
  ## Create mapping of values in additional columns to numeric values
  if (!is.null(additional_columns)) {
    additional_columns_legend <- vector(mode = "list", length = length(additional_columns)); names(additional_columns_legend) <- names(additional_columns)
    for (i in 1:length(additional_columns)) {
      name <- names(additional_columns)[i]
      vals <- unique(unlist(additional_columns[[i]]))
      flag <- if (is.numeric(vals)) { vals } else { 1:length(vals) }
      additional_columns_legend[[i]] <- data.frame(column = rep(name, length(vals)),
                                                   original = as.character(vals),
                                                   flag = flag)
    }
  }
  
  concat_efcs <- list()
  
  idcs.unused <- aberrances.list[[1]]$idcs.unused
  
  for (i in 1:length(fcs.paths)) {
    file      <- fcs.paths[i]
    file_used <- !(i %in% idcs.unused)
    
    cat("Enhancing ", basename(file), "\n", sep = "")
    if (!file_used) cat("(File was not used in analysis, but will be enhanced for easy concatenation of the whole dataset.)\n")
    
    ff <- read.FCS(file)
    
    if (is.null(junk)) {
      which.not_removed <- if (file_used) { 1:nrow(ff) } else { integer(0) }
    } else if (!file_used) {
      which.now_removed <- integer(0)
    } else {
      remove_idcs <- unique(unlist(junk[[basename(file)]]))
      
      which.not_removed <- (1:nrow(ff))
      if (length(remove_idcs) > 0) which.not_removed <- which.not_removed[-remove_idcs]
      
      nn1 <- length(which.not_removed)                           # clean events indices version 1 for checking
      a   <- aberrances.list[[1]]$binning[[basename(file)]]$assignment
      nn2 <- if (is.null(dim(a))) { length(a) } else { nrow(a) } # clean events indices version 2 for checking
      
      if (file_used && nn1 != nn2) stop("Indices misaligned")
    }
    
    ## Add event_used column
    event_used <- rep(0, nrow(ff))
    if (file_used) event_used[which.not_removed] <- 1000
    name       <- "event_used_in_analysis"
    cat("-> adding column ", name, "\n", sep = "")                          
    ff         <- fcs.add_col(ff, new_col = event_used,
                              colname = name)
    
    for (j in 1:length(aberrances.list)) {
      assignment          <- aberrances.list[[j]]$binning[[basename(file)]]$assignment
      if (!is.null(assignment) && is.null(dim(assignment))) assignment <- matrix(assignment)
      filter_responses    <- aberrances.list[[j]]$binning[[basename(file)]]$continuous
      
      if (!file_used) {
        if (aberrances.list[[j]]$method == "CellCnn") { # add dummy columns for concatenation
          nb               <- aberrances.list[[j]]$model$n_bins
          assignment       <- matrix(0, nrow = nrow(ff), ncol = nb)
          filter_responses <- matrix(0, nrow = nrow(ff), ncol = nb)
          # if (nb > 1) filter_responses_PC <- matrix(0, nrow = nrow(ff), ncol = nb)
        } else {
          assignment       <- matrix(0, nrow = nrow(ff), ncol = 1)
        }
      }
      
      ## Add binning & aberrance stats
      for (k in 1:ncol(assignment)) {
        name <- if (aberrances.list[[j]]$method == "CellCnn") { paste0(names(aberrances.list)[j], "_selection_", k, "_bin") } else { paste0(names(aberrances.list)[j], "_bin") }
        cat("-> adding column ", name, "\n", sep = "")
        binning <- rep(removed_flag, nrow(ff))
        if (file_used) binning[which.not_removed] <- assignment[, k]
        
        ff      <- fcs.add_col(ff, new_col = binning, colname = name)
        name <- if (aberrances.list[[j]]$method == "CellCnn") { paste0(names(aberrances.list)[j], "_selection_", k, "_aberrance") } else { paste0(names(aberrances.list)[j], "_aberrance") }
        ab_vals <- rep(0, length(binning))
        if (file_used) ab_vals[which.not_removed][assignment[, k] != 0] <- aberrances.list[[j]]$aberrance[[basename(file)]]$real[assignment[, k][assignment[, k] != 0]]
        
        overrepresentation                           <- ab_vals
        overrepresentation[overrepresentation < 0]   <- 0
        underrepresentation                          <- -ab_vals
        underrepresentation[underrepresentation < 0] <- 0
        
        m <- max(overrepresentation)
        if (m > 0) overrepresentation <- overrepresentation / m
        m <- max(underrepresentation)
        if (m > 0) underrepresentation <- underrepresentation / m
        
        name <- if (aberrances.list[[j]]$method == "CellCnn") { paste0(names(aberrances.list)[j], "_selection_", k, "_overrepresentation") } else { paste0(names(aberrances.list)[j], "_overrepresentation") }
        cat("-> adding column ", name, "\n", sep = "")
        ff <- fcs.add_col(ff, new_col = overrepresentation,
                          colname = name)
        name <- if (aberrances.list[[j]]$method == "CellCnn") { paste0(names(aberrances.list)[j], "_selection_", k, "_underrepresentation") } else { paste0(names(aberrances.list)[j], "_underrepresentation") }
        cat("-> adding column ", name, "\n", sep = "")
        ff <- fcs.add_col(ff, new_col = underrepresentation,
                          colname = name)
      }
      
      ## Add filter responses
      if (!is.null(filter_responses)) {
        for (k in 1:ncol(filter_responses)) {
          name <- paste0(names(aberrances.list)[j], "_filter_", k, "_response")
          cat("-> adding column ", name, "\n", sep = "")
          response <- rep(removed_flag, nrow(ff))
          if (file_used) response[which.not_removed] <- filter_responses[, k]
          
          ff <- fcs.add_col(ff, new_col = response, colname = name)
        }
      }
    }
    
    ## Dimension reduction
    if (!is.null(dimred.list) && length(dimred.list) > 0) {
      this_sample_idcs <- sample_idcs$clean[[i]][1]:sample_idcs$clean[[i]][2]
      cat("Sample idcs after filtering relevant events: ", this_sample_idcs[1], " ", this_sample_idcs[length(this_sample_idcs)], "\n", sep = "")
      
      for (dimred in dimred.list) {
        dimred_segment <- dimred$layout[this_sample_idcs, ]
        
        for (j in 1:ncol(dimred_segment)) {
          dimred_vec <- rep(removed_flag, nrow(ff))
          dimred_vec[which.not_removed] <- dimred_segment[, j]
          
          if (length(which.not_removed) != length(dimred_segment[, j])) stop("Dim reduction vector length error")
          
          name <- paste0(dimred$name, "_", j)
          cat("-> adding column ", name, "\n", sep = "")
          
          ff <- fcs.add_col(ff, new_col = dimred_vec, colname = name)
        }
      }
    }
    
    if (!file.exists(output.path)) {
      dir.create(output.path)
    }
    
    if (!is.null(additional_columns)) {
      for (j in 1:length(additional_columns)) {
        name  <- names(additional_columns)[j]
        value <- additional_columns[[j]][i]
        flag  <- additional_columns_legend[[j]]$flag[additional_columns_legend[[j]]$original == value]
        
        cat("-> adding custom column ", name, "\n", sep = "")
        
        ff    <- fcs.add_col(ff, new_col = rep(flag, nrow(ff)), colname = name)
      }
    }
    
    valid_ff <- make_valid_fcs(exprs = ff@exprs, desc1 = read.FCS(fcs.paths[1])@parameters$desc)
    write.FCS(valid_ff, filename = file.path(output.path, paste0(strsplit(basename(file), "[.]fcs")[[1]], "_enhanced.fcs")))
    
    if (!is.null(concatenated_fcs)) concat_efcs[[i]] <- cbind(ff@exprs, sample_idx = rep(i, nrow(ff@exprs)))
  }
  
  if (!is.null(additional_columns)) {
    additional_columns_legend <- do.call(rbind, additional_columns_legend)
    write.csv(additional_columns_legend, file.path(output.path, "_additional_columns_legend.csv"), row.names = FALSE)
  }
  
  sample_idcs_legend <- cbind(sample_idx = 1:length(fcs.paths), sample_name = basename(fcs.paths))
  write.csv(sample_idcs_legend, file.path(output.path, "_sample_idcs_legend.csv"), row.names = FALSE)
  
  if (!is.null(concatenated_fcs)) {
    concat_efcs <- do.call(rbind, concat_efcs)
    concat_ff   <- make_valid_fcs(exprs = concat_efcs, desc1 = read.FCS(fcs.paths[1])@parameters$desc)
    write.FCS(concat_ff, file.path(output.path, concatenated_fcs))
  }
}

