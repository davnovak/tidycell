## dimension reduction module

dimred.layout <- function(aa_s3,
                          channels_of_interest,
                          method.functor     = NULL,
                          method             = "tsne",
                          method.settings    = NULL,
                          normalise.quantile = NULL) {
    
    cols <- colnames(aa_s3$concat$exprs)
    efcs <- aa_s3$concat$exprs[, cols %in% channels_of_interest]
    
    if (!is.null(normalise.quantile)) {
        efcs <- apply(efcs, 2, function(col) {
            denom <- quantile(col, normalise.quantile)
            if (denom != 0) return(col/denom)
            col
        })
    }
    
    if (method == "tsne") {
        cat("Creating t-SNE projection...\n")
        initial_dims <- if (!is.null(method.settings$initial_dims)) { method.settings$initial_dims } else { 20 }
        num_threads  <- if (!is.null(method.settings$num_threads)) { method.settings$num_threads } else { 0 }
        
        dimred <- Rtsne::Rtsne(X                = efcs,
                               dims             = 2,
                               initial_dims     = initial_dims,
                               num_threads      = num_threads,
                               check_duplicates = FALSE)
        
        dimred$Y <- dimred$Y - min(dimred$Y) + 0.001
        
        return(list(name = "tsne",
                    layout = dimred$Y))
        
    } else if (method == "umap") {
        cat("Creating UMAP projection...\n")
        
        n_neighbours <- if (!is.null(method.settings$n_neighbours)) { method.settings$n_neighbours } else { 50L }
        
        dimred <- uwot::umap(efcs, n_neighbors = n_neighbours, n_components = 2)
        dimred <- dimred - min(dimred) + 0.001
        
        return(list(name = "umap",
                    layout = dimred))
        
    } else if (method == "embedsom") {
        cat("Creating EmbedSOM projection...\n")
        grid.x   <- if (!is.null(method.settings$grid.x))   { method.settings$grid.x }   else { 20 }
        grid.y   <- if (!is.null(method.settings$grid.y))   { method.settings$grid.y }   else { 20 }
        emcoords <- if (!is.null(method.settings$emcoords)) { method.settings$emcoords } else { "flat" }
        
        cat("emcoords =", emcoords, "\n")
        
        input     <- FlowSOM::ReadInput(flowFrame(efcs), scale = TRUE, silent = TRUE)
        som       <- FlowSOM::BuildSOM(input, xdim = grid.x, ydim = grid.y, colsToUse = 1:ncol(efcs), silent = TRUE)
        
        embedding <- suppressMessages(EmbedSOM::EmbedSOM(fsom = som, emcoords = emcoords))
        
        embedding <- embedding - min(embedding) + 0.001
        
        if (emcoords == "uwot::umap") emcoords <- "umap"
        
        return(list(name = paste0("embedsom_", emcoords),
                    layout = embedding))
    }
}
