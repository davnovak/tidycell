## stratified downsampling module


downsample <- function(labels,
                       N,
                       threshold    = N/100,
                       verbose      = TRUE,
                       ignore.label = -1,
                       seed         = NULL) { # returns indices in 'labels'
  
  ## To do: do this in C
  if (N > length(labels)) stop("Target size higher than original: cannot downsample")
  if(!is.null(seed)) set.seed(seed)
  
  labels.sorted <- as.numeric(names(sort(table(labels))))
  labels.sorted <- labels.sorted[labels.sorted != as.character(ignore.label)]
  
  if (N / length(labels.sorted) < threshold) warning("N too small for this threshold.")
  
  labels.inds <- lapply(labels.sorted, function(i) which(labels == i))
  
  out <- list()
  
  for (i in 1:length(labels.inds)) {
    if (verbose) cat(".")
    remaining_categories <- length(labels.inds) - i
    remaining_N          <- N - length(unlist(out))
    this_N               <- length(labels.inds[[i]])
    if (this_N < threshold) {
      if (verbose) cat("Population ", labels.sorted[i], " smaller than threshold: N = ", this_N, ".\n", sep = "")
      out[[i]] <- labels.inds[[i]]
    } else {
      selection_size <- min(remaining_N %/% (remaining_categories+1), this_N)
      out[[i]] <- sample(x = labels.inds[[i]], size = selection_size, replace = FALSE)
    }
  }
  
  if (verbose) cat("\n")
  out <- unlist(out)
  if (verbose) cat("Downsampled N = ", length(out), ".\n", sep = "")
  if (length(out) < N) warning("Downsampled N is lower than requested.")
  out
}
