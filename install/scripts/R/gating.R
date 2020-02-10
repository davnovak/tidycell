## gating module

## Function: get indices of a population from a vector of full paths or vector of short population names
.get_existing_pop_idcs <- function(query, pops, pops.short) {
  
  idx <- which(pops.short == query)
  if (length(idx) > 1) return(NULL)
  if (length(idx) == 0) {
    idx <- grep(paste0(query, "$"), pops)
    if (length(idx) > 1) return(NULL)
  }
  idx
}

## Function: get indices of parent's children populations from a vector of full paths of populations
.get_children_idcs <- function(query, pops) {
  
  all_idcs   <- grep(query, pops)
  parent_idx <- grep(paste0(query, "$"), pops[all_idcs])
  all_idcs[-parent_idx]
}

## Function: define a Boolean gate and retrieve indices of events within it.
## For instance, to get indices of all "Naive B cells" that are not "Transitional B cells", call
### get_new_pop_idcs("Naive B cells", exclude_children = "Transitional B cells")
## To get indices of events whose terminal node is "Viable cells", call
### get_new_pop_idcs("Viable cells", only_include_children = "")
## (That means exclude all children populations...)
## To get indices of events assigned to "B cells" or any of their subpopulations, simply call
### get_new_pop_idcs("B cells")
## In case any population is ambiguous (eg. suppose there are two nodes named "Terminal T cells")
## which are children of different parents -- "CD8 T cells" and "CD4 T cells"), then instead of the name
## of the population, include a substring of its path that unambiguous (eg. "CD8 T cells/Terminal T cells").
.get_new_pop_idcs <- function(parent,
                              exclude_children = NULL,
                              only_include_children = NULL,
                              pops = nodes,
                              idx_mat = idcs) {
  
  pops.short <- lapply(strsplit(pops, "/"), function(x) tail(x, 1))
  pops.short <- unlist(pops.short)
  
  parent_idx <- .get_existing_pop_idcs(parent, pops, pops.short)
  if (length(parent_idx) == 0) stop("Parent name non-existent or ambiguous")
  
  row_idcs   <- which(idx_mat[, parent_idx] > 0)
  
  if (!is.null(exclude_children)) {
    
    exclude_idcs <- lapply(exclude_children, function(x) .get_existing_pop_idcs(x, pops, pops.short))
    exclude_idcs <- unlist(exclude_idcs)
    
    submatrix    <- idx_mat[row_idcs, exclude_idcs, drop = FALSE]
    to_exclude   <- which(rowSums(submatrix) > 0)
    row_idcs     <- row_idcs[-to_exclude]
  }
  
  if (!is.null(only_include_children)) {
    
    remove_children <- (length(only_include_children) == 1 && only_include_children == "")
    
    if (remove_children)  {
      include_idcs <- .get_children_idcs(parent, pops)
    } else {
      include_idcs <- lapply(only_include_children, function(x) .get_existing_pop_idcs(x, pops, pops.short))
      include_idcs <- unlist(include_idcs)
    }
    
    submatrix  <- idx_mat[row_idcs, include_idcs, drop = FALSE]
    to_exclude <- if (remove_children) { which(rowSums(submatrix)) > 0 } else { which(rowSums(submatrix) == 0) }
    
    row_idcs <- row_idcs[-to_exclude]
  }
  
  row_idcs
}

## Function: get vector of per-event bin assignments for an .fcs file, provided the FlowJo workspace
gating_vector <- function(fcs.dir = "",
                          fcs.files,
                          wsp,
                          ungated_label = "*ungated*",
                          exclusive_pops = NULL,
                          new_Boolean_gates = NULL,
                          print_ids = TRUE) {
  
  #ws   <- CytoML::openWorkspace(wsp)
  ws   <- CytoML::open_flowjo_xml(wsp)
  #gs   <- CytoML::parseWorkspace(ws, name = 1, path = fcs.dir, sampNloc = "sampleNode")
  gs   <- flowjo_to_gatingset(ws, name = 1, path = fcs.dir, sampNloc = "sampleNode")
  #samp <- CytoML::getSampleGroups(ws)
  samp <- fj_ws_get_samples(ws)
  
  
  dd <- lapply(fcs.files, function(file) {
    ff <- read.FCS(file.path(fcs.dir, file))
    list(guid = ff@description$GUID,
         N    = nrow(ff@exprs))
  })
  
  ids <- list()
  for (name in sampleNames(gs)) {
    for (i in 1:length(gs)) {
      if (gs[[i]]@name == name) {
        ids[[i]] <- name
        break()
      }
    }
  }
  ids <- unlist(ids)
  
  out <- vector(mode = "list", length = length(ids)); names(out) <- fcs.files
  for (i in 1:length(out)) {
    id    <- ids[i]
    
    #nodes <- flowWorkspace::getNodes(gs[[id]])[-1]
    nodes <- CytoML::gs_get_pop_paths(gs[[id]])[-1]
    #idcs  <- flowWorkspce::getIndiceMat(gs[[id]], paste(nodes, collapse = "|")) # full indices matrix
    idcs  <- CytoML::gh_pop_get_indices_mat(gs[[id]], paste(nodes, collapse = "|")) # full indices matrix
    
    N     <- nrow(idcs)
    
    if (!is.null(new_Boolean_gates)) { # new_Boolean_gates overrides exclusive_pops
      gates <- rep(ungated_label, N)
      
      for (j in 1:length(new_Boolean_gates)) {
        name                  <- names(new_Boolean_gates)[j]
        new_gate              <- new_Boolean_gates[[j]]
        parent                <- new_gate$parent
        exclude_children      <- new_gate$exclude_children
        only_include_children <- new_gate$only_include_children
        
        gates[.get_new_pop_idcs(parent                = parent,
                                exclude_children      = exclude_children,
                                only_include_children = only_include_children,
                                pops                  = nodes,
                                idx_mat               = idcs)] <- name
      }
      
      out[[i]] <- gates
      
    } else {
      if (!is.null(exclusive_pops)) {
        relevant <- lapply(exclusive_pops, function(x) grep(paste0("/", x), nodes))
        
        idcs.contracted <- lapply(relevant, function(x) { # contracted indices matrix
          i                   <- idcs[, x, drop = FALSE]
          vec                 <- rep(0, N)
          vec[rowSums(i) > 0] <- 1
          vec
        }) %>% do.call(cbind, .)
        
        colnames(idcs.contracted) <- exclusive_pops
        
        idcs <- idcs.contracted
      }
      
      gates                            <- colnames(idcs)[apply(idcs, MARGIN = 1, which.max)]
      gates[which(rowSums(idcs) == 0)] <- ungated_label
      out[[i]]                         <- gates
    }
  }
  
  out
}

## Function: define custom populations based on a combination of positive & negative phenotypes in specified
## parameters. Get sizes of these populations in controls & affected, compare and report Wilcoxon test p-value.
custom_population_counts <- function(fcs.paths,
                                     fcs.files,
                                     idcs.controls,
                                     idcs.affected,
                                     junk,
                                     thresholds,
                                     use_marker_names = TRUE,
                                     channels.positive = NULL,
                                     channels.negative = NULL) {
  
  .get_pop_counts <- function(idx, path) {
    ff   <- read.clean_FCS(path, junk, idx)
    efcs <- ff@exprs
    if (use_marker_names) colnames(efcs) <- FlowSOM::get_markers(ff, markers = colnames(ff))
    rm(ff)
    
    positive <- if (is.null(channels.positive)) { rep(TRUE, nrow(efcs)) } else {
      Reduce("&", lapply(channels.positive, function(y) efcs[, y] > thresholds[[idx]][y]))
    }
    
    negative <- if (is.null(channels.negative)) { rep(TRUE, nrow(efcs)) } else {
      Reduce("&", lapply(channels.negative, function(y) efcs[, y] < thresholds[[idx]][y]))
    }
    
    sum(positive & negative)/nrow(efcs)
  }
  
  vals0        <- sapply(1:length(fcs.paths[idcs.controls]), function(i) .get_pop_counts(idx = i, path = fcs.paths[idcs.controls][i]))
  names(vals0) <- basename(fcs.paths[idcs.controls])
  
  vals1        <- sapply(1:length(fcs.paths[idcs.affected]), function(i) .get_pop_counts(idx = i, path = fcs.paths[idcs.affected][i]))
  names(vals1) <- basename(fcs.paths[idcs.affected])
  
  suppressWarnings(pv <- wilcox.test(vals0, vals1, alternative = "two.sided")$p.value)

  list(counts.controls = vals0,
       counts.affected = vals1,
       p.value         = pv)
}

get_gate_indices <- function(fcs.dir,
                             fcs.files,
                             wsp,
                             gate_names,
                             invert = FALSE) {
  
  #ws   <- CytoML::openWorkspace(wsp)
  ws   <- CytoML::open_flowjo_xml(wsp)
  #gs   <- CytoML::parseWorkspace(ws, name = 1, path = fcs.dir, sampNloc = "sampleNode")
  gs   <- CytoML::flowjo_to_gatingset(ws, name = 1, path = fcs.dir)
  #samp <- CytoML::getSampleGroups(ws)
  samp <- CytoML::fj_ws_get_samples(ws)
  
  
  dd <- lapply(fcs.files, function(file) {
    ff <- read.FCS(file.path(fcs.dir, file))
    list(guid = ff@description$GUID,
         N    = nrow(ff@exprs))
  })
  
  ids <- list()
  for (name in sampleNames(gs)) {
    for (i in 1:length(gs)) {
      if (gs[[i]]@name == name) {
        ids[[i]] <- name
        break()
      }
    }
  }
  ids <- unlist(ids)
  
  out <- vector(mode = "list", length = length(ids)); names(out) <- fcs.files
  for (i in 1:length(out)) {
    
    id <- which(gsub("[.]fcs_[0-9]+", ".fcs", flowWorkspace::sampleNames(gs)) == names(out)[i])
    
    cat("Sample ID:", id, "\n")
    cat("Match:", samp[samp$sampleID == id, ]$name, gs[[i]]@name, names(out)[i], sep = "\n")
    
    #nodes <- flowWorkspace::getNodes(gs[[id]])[-1]
    nodes <- flowWorkspace::gs_get_pop_paths(gs[[id]])[-1]
    
    #idcs  <- flowWorkspace::getIndiceMat(gs[[id]], paste(nodes, collapse = "|"))
    idcs  <- flowWorkspace::gh_pop_get_indices_mat(gs[[id]], paste(nodes, collapse = "|"))
    
    cat("N_events in .wsp: ", nrow(idcs), "\nN_events in .fcs: ", dd[[i]]$N, "\n\n", sep = "")
    
    which.cols <- unique(unlist(lapply(gate_names, function(gate_name) which(stringr::str_detect(colnames(idcs), stringr::fixed(gate_name))))))
    
    idcs         <- as.matrix(idcs[, which.cols], ncol = 1)
    which.events <- which(rowSums(idcs) > 0)
    
    if (invert) which.events <- if (length(which.events) > 0) { (1:nrow(idcs))[-which.events] } else { 1:nrow(idcs) }
    
    cat("Selected ", length(which.events), " out of ", nrow(idcs), " events (complement = ", nrow(idcs) - length(which.events),").\n\n", sep = "")
    
    out[[i]] <- which.events
  } 
  
  out
}
