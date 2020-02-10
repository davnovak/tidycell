## input check module

input_check <- function(csv.analysis_inputs = "analysis_inputs.csv",
                        csv.analysis_markers = "analysis_markers.csv",
                        fcs.directory,
                        file_connection = NULL) {
  
  msg        <- vector(mode = "character")
  fcs_files  <- as.character(read.table(csv.analysis_inputs, sep = ",", header = TRUE)$fcs_file)
  
  m          <- read.table(csv.analysis_markers, sep = ",", header = TRUE)
  params.csv <- data.frame(channel = as.character(m$channels),
                           marker = as.character(m$markers),
                           stringsAsFactors = FALSE)
  
  fcs_paths  <- file.path(fcs.directory, fcs_files)
  
  errors     <- FALSE
  
  for (path in fcs_paths) {
    str <- paste("Checking file: ", basename(path), "...\n", sep = ""); cat(str)
    msg <- c(msg, str)
    
    if (!file.exists(path)) {
      str    <- paste("-> file not found\n", sep = ""); cat(str, "\n")
      msg    <- c(msg, str)
      errors <- TRUE
    } else {
      ff     <- read.FCS(path)
      desc   <- na.omit(ff@parameters$desc)
      
      channels <- sapply(names(desc), function(marker_slot) {
        channel_slot <- gsub(".{1}$", "N", marker_slot)
        channel      <- ff@description[[channel_slot]]
      })
      markers <- as.vector(desc)
      
      params.ff <- data.frame(channel = channels,
                              marker = markers,
                              stringsAsFactors = FALSE)
      
      missing_markers  <- params.csv$marker[sapply(params.csv$marker, function(m) !m %in% params.ff$marker)]
      missing_channels <- params.csv$channel[sapply(params.csv$channel, function(m) !m %in% params.ff$channel)]
      if (length(missing_markers) > 0) {
        str    <- paste("-> missing markers: ", paste(missing_markers, collapse = ", ")); cat(str, "\n")
        msg    <- c(msg, str)
        errors <- TRUE
      }
      if (length(missing_channels) > 0) {
        str    <- paste("-> missing channels: ", paste(missing_channels, collapse = ", ")); cat(str, "\n")
        msg    <- c(msg, str)
        errors <- TRUE
      }
      
      for (idx.channel in 1:length(params.csv$channel)) {
        channel.csv <- params.csv$channel[idx.channel]
        marker.csv  <- params.csv$marker[idx.channel]
        
        marker.ff   <- params.ff$marker[params.ff$channel == channel.csv]
        channel.ff  <- channel.csv
        
        if (length(marker.ff) == 0 || length(channel.ff) == 0) next()
        
        if (marker.csv != marker.ff) {
          str    <- paste("-> mismatch: ", channel.ff, " paired with ", marker.ff, ", expected ", marker.csv, sep = ""); cat(str, "\n")
          msg    <- c(msg, str)
          errors <- TRUE
        }
      }
    }
  }
  
  msg <- c(msg, "\n")
  
  if (!is.null(file_connection)) writeLines(msg, file_connection)
  
  cat("Done\n")
  errors
}