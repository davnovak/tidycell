## Extract channels and markers for analysis

extract_markers_input <- function(fcs_directory,
                                  output_name = "markers_and_channels.csv") {
  
  if (file.exists(output_name)) stop("File with output name already exists! Please choose a different name")
  
  files    <- list.files(fcs_directory, pattern = ".fcs", full.names = TRUE)
  
  channels <- list()
  markers  <- list()
  
  for (i in 1:length(files)) {
    file <- files[i]
    ff   <- flowCore::read.FCS(file)
    
    desc <- na.omit(ff@parameters$desc)
    
    channels[[i]] <- sapply(names(desc), function(marker_slot) {
      channel_slot <- gsub(".{1}$", "N", marker_slot)
      channel      <- ff@description[[channel_slot]]
    })
    markers[[i]] <- as.vector(desc)
  }
  
  names(channels) <- names(markers) <- n <- basename(files)
  
  out <- lapply(n, function(x) {
    data.frame(channels = channels[[x]],
               markers = markers[[x]])
  })
  names(out) <- n
  out        <- do.call(cbind, out)
  
  write.table(out, file = output_name, sep = ",", row.names = FALSE)
}

## USAGE EXAMPLE:
# extract_markers_input("fcs_files/", "my_channels.csv")
