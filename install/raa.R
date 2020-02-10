package_list <- c('flowCore',
                  'flowWorkspace',
                  'tidyverse',
                  'Rtsne',
                  'uwot',
                  'FlowSOM',
                  'flowDensity',
                  'flowAI',
                  'CytoML',
                  'RColorBrewer',
                  'cowplot',
                  'pheatmap',
                  'gplots',
                  'scales',
                  'reticulate')
for (pkg in package_list) suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))

use_condaenv("tidycell_cellcnn_env")

tidycell_path <- "/home/rstudio/tidycell/scripts"

for (file in list.files(file.path(tidycell_path, "R"), full.names = TRUE)) source(file)

setwd("./data")

cat("#### tidycell aberrance analysis ####\n")

call.aberrance_analysis(csv.analysis_inputs   = "analysis_inputs.csv",
                        csv.analysis_markers  = "analysis_markers.csv",
                        csv.analysis_settings = "analysis_settings.csv",
                        tidycell_path = tidycell_path)
