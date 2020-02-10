## reindexing module

.reindex <- function(controls, affected, unused) {
    
    if (is.null(unused)) {
        return(list(controls = controls,
                    affected = affected))
    }
    concat <- c(controls, affected, unused)
    labels <- c(rep(0, length(controls)),
                rep(1, length(affected)),
                rep(-1, length(unused)))
    
    labels <- labels[order(concat)]
    labels <- labels[labels != -1]
    
    list(controls = which(labels == 0),
         affected = which(labels == 1))
}