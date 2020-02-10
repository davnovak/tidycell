## channel extraction module

## Function: extract channel and marker names from a flowFrame, with the ability to exclude some channels
## by the name of channel, name of marker or column index.
## (Does not include scatter or time parameters.)
fcs_params <- function(ff,
                       exclude_channels = NULL,
                       exclude_markers  = NULL,
                       exclude_idcs     = NULL) {
  
  markers_of_interest  <- flowCore::markernames(ff)
  channels_of_interest <- colnames(ff)
  
  channels_of_interest       <- channels_of_interest[!channels_of_interest %in% c("FSC-A", "FSC-W", "FSC-H", "SSC-A", "SSC-W", "SSC-H", "Time")]
  names(markers_of_interest) <- channels_of_interest
  
  idcs         <- 1:length(markers_of_interest)
  exclude_idcs <- unique(c(exclude_idcs, which(markers_of_interest %in% exclude_markers, which(channels_of_interest %in% exclude_channels))))
  
  if (length(exclude_idcs) > 0) {
    channels_of_interest <- channels_of_interest[-exclude_idcs]
    markers_of_interest  <- markers_of_interest[-exclude_idcs]
    idcs                 <- idcs[-exclude_idcs]
  }
  
  list(channels = channels_of_interest,
       markers = markers_of_interest,
       idcs = idcs)
}