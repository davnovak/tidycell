## visualisation module

## Function: create plot showing composition composition of bins in terms of gated populations
create_visual.bin_composition <- function(aberrances,
                                          sample_idx,
                                          gates,
                                          only_show_significant = FALSE,
                                          alpha = 0.1,
                                          visual.path = NULL) {
  
  n_bins  <- aberrances$n_bins
  gates   <- unlist(gates)
  gates   <- clean_vector(gates, junk, sample_idx)
  binning <- aberrances$binning[[sample_idx]]$assignment
  gates.u <- unique(gates)
  
  which_bins <- if (only_show_significant) { which(aberrances$p.vals < alpha) } else { 1:n_bins }
  
  gate_representation <- lapply(which_bins, function(bin) {
    gates_in_bin           <- data.frame(table(gates[binning == bin]))
    missing_gates          <- gates.u[!gates.u %in% gates_in_bin[, 1]]
    gates_in_bin           <- rbind(gates_in_bin, data.frame(Var1 = missing_gates,
                                                             Freq = rep(0, length(missing_gates))))
    gates_in_bin           <- cbind(rep(bin, nrow(gates_in_bin)), gates_in_bin)
    colnames(gates_in_bin) <- c("bin", "gate", "count")
    
    gates_in_bin[, 3]      <- gates_in_bin[, 3] / sum(gates_in_bin[, 3])
    gates_in_bin
  })
  gate_representation      <- do.call(rbind, gate_representation)
  gate_representation$bin  <- as.factor(gate_representation$bin)
  gate_representation$gate <- as.factor(gate_representation$gate)
  
  suppressWarnings(
    palette_getter <- colorRampPalette(brewer.pal(length(gates.u), "Set1"))
  )
  colours          <- palette_getter(length(gates.u))
  colours[1]       <- "grey" # corresponds to ungated events
  
  bin_composition_plot <- ggplot(data = gate_representation,
                                 aes(x = bin, y = count)) +
    geom_col(aes(fill = gate), width = 0.8) +
    ylab("proportion in bin") +
    theme(panel.background = element_blank(), axis.line = element_blank()) +
    scale_fill_manual(values = colours) +
    ggtitle("Bin composition in terms of gated populations and plurality populations per bin")
  
  largest_populations_labels <- lapply(which_bins, function(bin) {
    representation  <- gate_representation[gate_representation$bin == bin, ]
    which.plurality <- which.max(representation$count)
    plurality.gate  <- representation$gate[which.plurality]
    plurality.count <- round(representation$count[which.plurality] * 100, 2)
    return(paste0("Bin ", bin, ": ", plurality.gate, " (", plurality.count, "%)"))
  })
  
  largest_populations_labels <- paste0(largest_populations_labels, collapse = "\n")
  
  largest_populations_plot   <- ggplot() + annotate("text", x = 4, y = 25, size = 4, label = largest_populations_labels) + theme_bw() +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.line = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), legend.position = "none",
          panel.background = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  bin_composition_grid <- plot_grid(bin_composition_plot, largest_populations_plot, align = "h", ncol = 2, rel_widths = c(5/7, 2/7))
  
  if (!is.null(visual.path)) {
    png(visual.path,
        width  = 1200,
        height = 900)
    plot(bin_composition_grid)
    dev.off()
  }
  
  bin_composition_grid
}

## Function: provided a FlowSOM binning model, plot pie-charts showing expressions of specified marker combinations
create_visual.som_stars <- function(aberrances,
                                    title,
                                    markers,
                                    visual.path = NULL) {
  
  if (aberrances$method != "FlowSOM") stop("Cannot use PlotStars for other than FlowSOM model")
  
  if (!is.null(visual.path)) png(visual.path, width = 900, height = 900)
  
  PlotStars(aberrances$model,
            backgroundValues = aberrances$model$metaclustering,
            markers          = which(aberrances$model$prettyColnames %in% markers),
            main             = title,
            view             = "grid")
  
  if (!is.null(visual.path)) dev.off()
}

## Function: plot numbers of events falling into different bins and the difference therein between controls and affected
create_visual.bin_sizes <- function(aberrances,
                                    only_show_significant = FALSE,
                                    alpha = 0.1,
                                    bin.order = NULL,
                                    title,
                                    visual.path = NULL) {
  
  n_bins     <- aberrances$n_bins
  which_bins <- 1:n_bins
  if (only_show_significant) which_bins <- which(aberrances$p.vals < alpha)
  
  total.controls <- lapply(which_bins, function(bin) {
    lapply(aberrances$idcs.controls, function(idx) {
      assignment <- aberrances$binning[[idx]]$assignment
      sum(assignment == bin)/length(assignment)
    }) %>% unlist %>% as.numeric
  })
  
  total.affected <- lapply(which_bins, function(bin) {
    lapply(aberrances$idcs.affected, function(idx) {
      assignment <- aberrances$binning[[idx]]$assignment
      sum(assignment == bin)/length(assignment)
    }) %>% unlist %>% as.numeric
  })
  
  median.controls <- sapply(total.controls, median) * 100 %>% round(2)
  median.affected <- sapply(total.affected, median) * 100 %>% round(2)
  mad.controls    <- sapply(total.controls, mad) * 100 %>% round(2)
  mad.affected    <- sapply(total.affected, mad) * 100 %>% round(2)
  
  size_changes.median <- data.frame(bin        = as.factor(rep(which_bins, 2)),
                                    category   = c(rep("controls", length(which_bins)), rep("affected", length(which_bins))),
                                    percentage = c(median.controls, median.affected))
  size_changes.mad <- data.frame(bin        = as.factor(rep(which_bins, 2)),
                                 category   = c(rep("controls", length(which_bins)), rep("affected", length(which_bins))),
                                 percentage = c(mad.controls, mad.affected))
  
  bar_plot <- ggplot(data = size_changes.median,
                     aes(x = bin, y = percentage, group = factor(category, levels = c("controls", "affected")))) +
    geom_bar(stat = "identity", aes(fill = category), position = "dodge") +
    geom_errorbar(stat = "identity",
                  aes(ymin = size_changes.median$percentage - size_changes.mad$percentage,
                      ymax = size_changes.median$percentage + size_changes.mad$percentage),
                  colour = "darkblue",
                  size = 0.5,
                  position = "dodge") +
    labs(title = title,
         subtitle = "Controls versus affected (median and MAD are shown)") +
    ylab("percentage") +
    theme_minimal()
  
  if (!is.null(visual.path)) {
    png(visual.path,
        width  = 1200,
        height = 900)
    plot(bar_plot)
    dev.off()
  }
  
  bar_plot
}

## Function: plot bin aberrance (measure of how under- or over-represented each bin is) in controls versus affected
create_visual.aberrance <- function(aberrances,
                                    only_show_significant = FALSE,
                                    alpha = 0.1,
                                    bin.order = NULL,
                                    title,
                                    average.measure = "median",
                                    visual.path = NULL) {
  
  idcs.controls <- aberrances$idcs.controls
  idcs.affected <- aberrances$idcs.affected
  
  for_plot.controls <- lapply(aberrances$aberrance[idcs.controls], function(x) x$plot)
  for_plot.affected <- lapply(aberrances$aberrance[idcs.affected], function(x) x$plot)
  
  aberrance.controls <- lapply(aberrances$aberrance[idcs.controls], function(x) x$real)
  aberrance.affected <- lapply(aberrances$aberrance[idcs.affected], function(x) x$real)
  
  d0      <- lapply(aberrance.controls, function(x) data.frame(bin = as.factor(1:length(x)), aberrance = x))
  d0.plot <- lapply(for_plot.controls,  function(x) data.frame(bin = as.factor(1:length(x)), aberrance = x))
  d1      <- lapply(aberrance.affected, function(x) data.frame(bin = as.factor(1:length(x)), aberrance = x))
  d1.plot <- lapply(for_plot.affected,  function(x) data.frame(bin = as.factor(1:length(x)), aberrance = x))
  d0      <- do.call(rbind, d0)
  d0.plot <- do.call(rbind, d0.plot)
  d1      <- do.call(rbind, d1)
  d1.plot <- do.call(rbind, d1.plot)
  p.vals  <- aberrances$p.vals
  d0      <- d0.plot
  d1      <- d1.plot
  rm(d0.plot)
  rm(d1.plot)
  
  stats0 <- lapply(sort(unique(d0$bin)), function(bin) {
    vals <- d0$aberrance[d0$bin == bin]
    list(average   = median(vals),
         deviation = mean(vals))
  })
  stats0 <- data.frame(cbind(bin       = 1:length(stats0),
                             average   = sapply(stats0, function(x) x$average),
                             deviation = sapply(stats0, function(x) x$deviation)))
  
  stats1 <- lapply(sort(unique(d1$bin)), function(bin) {
    vals <- d1$aberrance[d1$bin == bin]
    list(average = median(vals),
         deviation = mad(vals))
  })
  stats1 <- data.frame(cbind(bin       = 1:length(stats1),
                             average   = sapply(stats1, function(x) x$average),
                             deviation = sapply(stats1, function(x) x$deviation)))
  
  d0 <- cbind(d0, category = rep("control", nrow(d0)))
  d1 <- cbind(d1, category = rep("affected", nrow(d1)))
  d  <- rbind(d0, d1)
  
  if (is.null(bin.order)) bin.order <- 1:length(aberrance.affected[[1]])
  
  which.significant <- which(p.vals < alpha)
  s                 <- ifelse(d$bin %in% which.significant, 1, 0)
  s                 <- as.factor(s)
  levels(s)         <- c("0", "1")
  d                 <- cbind(d, significance = s)
  
  if (only_show_significant) d <- d[d$bin %in% which.significant, ]
  
  p <- ggplot(d, aes(x = factor(bin, levels = bin.order), y = aberrance)) +
    geom_point(aes(colour = category), position = position_jitterdodge(), size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgreen") +
    theme_minimal() + ggtitle(title) +
    scale_color_manual(values = c("#0868ac", "#e34a33")) +
    xlab("bin")
  
  if (!only_show_significant) {
    p <- p + scale_fill_manual(values = c("white", "#fff7bc")) +
      geom_boxplot(aes(colour = category, fill = significance), position = position_dodge2(), outlier.shape = NA)
  } else {
    p <- p + geom_boxplot(aes(colour = category), position = position_dodge2(), outlier.shape = NA)
  }
  
  if (!is.null(visual.path)) {
    png(visual.path,
        width = 1200,
        height = 900)
    plot(p)
    dev.off()
  }
  
  p
}

## Function: create a collection of 2D scatterplots showing phenotype of cells in a selected bin
create_visual.bin_scatterplots <- function(aggregate,
                                           aberrances,
                                           which.bin,
                                           markers = NULL,
                                           markers_of_interest = NULL,
                                           alternative_colnames = NULL,
                                           visual.path.png = paste0("bin", which.bin, "_scatterplots.png"),
                                           visual.path.pdf = paste0("bin", which.bin, "_scatterplots.pdf")) {
  
  if (is.null(markers)) {
    if (is.null(markers_of_interest)) {
      stop("Either put 'markers' as list of pairwise combinations of markers or provide 'markers_of_interest' if you want all pairwise combinations to be plotted")
    }
    combos  <- combn(markers_of_interest, 2)
    markers <- split(combos, rep(1:ncol(combos), each = nrow(combos)))
  }
  
  expression_matrix <- aggregate$aggregate
  if (!is.null(alternative_colnames)) colnames(expression_matrix) <- alternative_colnames
  
  downsample_idcs   <- sample(1:nrow(expression_matrix), 10000)
  binning           <- aberrances$aggregate_matrix_binning
  idcs              <- if (!aberrances$bin_overlap) { which(binning == which.bin) } else { which(binning[, which.bin] > 0) }
  idcs              <- sample(idcs, min(length(idcs), 1000))
  
  plot_fun <- function() {
    for (marker in markers) {
      plot(expression_matrix[downsample_idcs, marker], pch = 20, cex = 4, col = "grey")
      points(expression_matrix[idcs, marker], pch = 20, cex = 4, col = rgb(red = 0.5, green = 0, blue = 0, alpha = 0.4))
    }
  }
  
  if (!is.null(visual.path.png)) {
    png(visual.path.png,
        width = 2 * 500,
        height = ceiling(length(markers) / 2) * 500)
    layout(matrix(1:(ceiling(length(markers) / 2) * 2), ncol = 2, byrow = TRUE))
    plot_fun()
    dev.off()
  }
  
  if (!is.null(visual.path.pdf)) {
    pdf(visual.path.pdf)
    plot_fun()
    dev.off()
  }
}

create_visual.bin_dimred <- function(aberrances,
                                     dimred_layout,
                                     sample_idx,
                                     bins = NULL,
                                     sleep = 0,
                                     plot.title = "",
                                     visual.path = NULL) {
  
  assignment <- aberrances$binning[[sample_idx]]$assignment
  is_mat     <- !is.null(dim(assignment))
  N          <- ifelse(!is_mat, length(assignment), nrow(assignment))
  all_bins   <- 1:max(assignment)
  if (is.null(bins)) bins <- all_bins
  if (max(bins) > 1) palette(rainbow(max(bins)))
  
  if (!is.null(visual.path)) {
    png(visual.path,
        width  = 900,
        height = 900)
  } else {
    X11()
  }
  
  plot(dimred_layout, col = "grey", pch = 20, axes = FALSE, xlab = "", ylab = "", cex = 1,
       main = plot.title)
  for (bin in bins) {
    Sys.sleep(sleep)
    
    if (!bin %in% all_bins) {
      warning(paste0("Bin ", bin, " not found"))
      continue()
    }
    idcs <- if (is_mat) { which(assignment[, bin] > 0) } else { which(assignment == bin) }
    points(dimred_layout[idcs, ], col = alpha(bin + 1, .4), pch = 20)
  }
  
  if (!is.null(visual.path)) dev.off()
}

## Function: plot a heatmap showing the expression profile of each bin (in terms of median fluorescence intensities)
create_visual.bin_heatmap <- function(aberrances,
                                      only_show_significant = FALSE,
                                      alpha = 0.1,
                                      out.bin_order = NULL,
                                      out.marker_order = NULL,
                                      visual.path = NULL,
                                      scale = "column") {
  
  MFIs   <- aberrances$MFIs
  if (only_show_significant) MFIs <- MFIs[, which(aberrances$p.vals < alpha)]
  visual <- pheatmap(MFIs, scale = scale)
  if (!is.null(visual.path)) {
    png(filename = visual.path,
        width    = 900,
        height   = 900)
    print(visual)
    dev.off()
  }
  if (!is.null(out.bin_order))    out.bin_order    <- visual$tree_col[["order"]]
  if (!is.null(out.marker_order)) out.marker_order <- visual$tree_row[["order"]]
  visual
}

create_visual.mem_scores <- function(scores,
                                     idcs.controls,
                                     idcs.affected,
                                     idcs.unused = NULL,
                                     visual.path = NULL) {
  require(ggplot2)
  
  if (!is.null(idcs.unused)) {
    relative_idcs <- .reindex(idcs.controls, idcs.affected, idcs.unused)
    idcs.controls <- relative_idcs$controls
    idcs.affected <- relative_idcs$affected
  }
  
  plot.per_sample <- lapply(1:length(scores$scores.per_sample), function(idx) {
    category <- ifelse(idx %in% idcs.controls, "controls", "affected")
    vals     <- scores$scores.per_sample[[idx]]
    data.frame(category = rep(category, length(vals)),
               marker   = names(vals),
               score    = vals)
  })
  plot.per_sample          <- do.call(rbind, plot.per_sample)
  plot.per_sample$category <- factor(plot.per_sample$category)
  
  visual <- ggplot(data = plot.per_sample, aes(x = marker, y = score)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
    geom_boxplot(aes(colour = category), position = position_dodge2(), outlier.shape = NA) +        
    geom_point(aes(colour = category), position = position_jitterdodge(), size = 0.6) +
    ggtitle("MEM scores for each marker") +
    theme_minimal()
  
  if (!is.null(visual.path)) {
    png(visual.path,
        width = 1200,
        height = 900)
    plot(visual)
    dev.off()
  }
  
  visual
}
