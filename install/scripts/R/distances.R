## distance function module

.dist.correlation <- function(x, y) cor(x, y, method = "pearson")

distance_matrix <- function(df, dist_fun = "correlation") {
  
  n <- ncol(df)
  
  if (dist_fun == "correlation") {
    .dist <- .dist.correlation
  } else {
    stop("Invalid distance function")
  }
  
  d <- matrix(nrow = n, ncol = n, dimnames = list(colnames(df), colnames(df)))
  for (i in 1:n)
    for (j in 1:i)
      if (i != j)
        d[i, j] <- .dist(df[, i], df[, j])
  
  as.dist(d)
}
