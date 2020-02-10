aa <- function() {
    aa_env             <- new.env(hash = TRUE)
    aa_env$concat      <- NULL
    aa_env$agg         <- NULL
    aa_env$dimred.list <- NULL
    structure(aa_env, class = 'aberrance_analysis_data')
}