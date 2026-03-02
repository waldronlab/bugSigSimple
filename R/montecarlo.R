#' Identify the most frequently recurring taxa in bugsigdb
#'
#' @param sigs data.frame produced by \link[bugsigdbr]{importBugSigDB}
#' @param n Number of top taxa to return
#'
#' @return integer vector of top most recurrent taxa
#' @export
#'
#' @examples
#' dat <- bugsigdbr::importBugSigDB()
#' dat.select <- bugsigdbr::getSignatures(dat, tax.level = "genus", exact.tax.level=TRUE)
#' dat.select <- dat.select[grep("UP", names(dat.select))] #only "UP" signatures
#' frequencySigs(sigs=dat.select, n = 10)
frequencySigs <- function(sigs, n = 10){
  sigs.tab <- sort(table(unlist(sigs)), decreasing=TRUE)
  names(sigs.tab) <- vapply(names(sigs.tab), .getTip,
                            character(1), USE.NAMES = FALSE)
  head(sigs.tab, n=n) 
}


#' Simulate a list of signatures based on a universe of taxa equal in size and length to some smaller set of signatures
#'
#' @param relevant.sigs A list of signatures providing a relevant background for drawing taxa
#' @param siglengths An integer vector, the length of which provides the number of signatures to be simulated, and the 
#' integers of which provide the number of taxa to be simulated in each signature
#'
#' @return A list of signatures of the same number and individual lengths as found in my.dat
#' @export
#'
#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' my.dat <- full.dat[full.dat$Curator == "Mst Afroza Parvin", ]
#' relevant.dat <- full.dat[full.dat$`Body site` %in% my.dat$`Body site`, ]
#' relevant.sigs <- bugsigdbr::getSignatures(relevant.dat, tax.level = "genus")
#' siglengths <- sapply(my.dat, length)
#' simulateSignatures(relevant.sigs, siglengths)

#simulateSignatures <-
#  function(relevant.sigs,
  #         siglengths) {
  #  sigs.universe <- unlist(relevant.sigs)
  #  universetable <- table(sigs.universe)
   # universetable <- universetable / sum(universetable)
  #  lapply(siglengths, function(n)
   #   sample(names(universetable), size = n, prob = universetable))
  #}

#.countBug <- function(relevant.sigs, siglengths){
#  siglist <- simulateSignatures(relevant.sigs, siglengths)
 # max(table(unlist(siglist)))
#}
simulateSignatures <- function(relevant.sigs, siglengths) {
  sigs.universe  <- unlist(relevant.sigs)
  universetable  <- table(sigs.universe)
  universetable  <- universetable / sum(universetable)
  lapply(siglengths, function(n)
    sample(names(universetable), size = n, prob = universetable))
}

.countBug <- function(relevant.sigs, siglengths) {
  siglist <- simulateSignatures(relevant.sigs, siglengths)
  max(table(unlist(siglist)))
}

#' countBug counts the frequency of the most commonly identified bug in a simulated signature.

#' getCriticalN performs a Monte Carlo simulation to estimate the number of times the most frequent taxon is expected to be observed
#' in a list of signatures
#'
#' @param relevant.sigs a list of signatures that form the "background" from which taxa for simulated signatures will be drawn. 
#' These are used to estimate how frequently taxa occur
#' @param siglengths The sizes of signatures found in a set of related studies. Simulated signatures will match these in number and size.
#' @param alpha Probability at which a critical threshold will be calculated (default: 0.05)
#' @param nsim Number of simulations (default: 1000)
#'
#' @return The 1 - alpha quantile of Monte Carlo simulated values for the maximum number of times any taxon is identified.
#' @export
#' @details E.g. for alpha = 0.05, we expect only a 5% chance that any taxon will be identified N times or more.

#' @examples
#' full.dat <- bugsigdbr::importBugSigDB()
#' my.dat <- full.dat[full.dat$Curator == "Mst Afroza Parvin", ]
#' relevant.dat <- full.dat[full.dat$`Body site` %in% my.dat$`Body site`, ]
#' relevant.sigs <- bugsigdbr::getSignatures(my.dat)
#' my.sigs.increased <- relevant.sigs[grep("UP", names(relevant.sigs))]
#' (my.siglengths <- sapply(my.sigs.increased, length))
#' getCriticalN(relevant.sigs, my.siglengths)
#' # Compare to observed
#' frequencySigs(my.sigs.increased)

#getCriticalN <- function(relevant.sigs, siglengths, alpha = 0.05, nsim = 1000){
#  res <- replicate(nsim, suppressMessages(.countBug(relevant.sigs, siglengths)))
 # quantile(res, 1 - alpha)
#}
getCriticalN <- function(relevant.sigs,
                         siglengths,
                         alpha    = 0.05,
                         nsim     = 1000,
                         ci       = 0.95,
                         ci_nsim  = 500) {
  
  sim_vals   <- replicate(nsim,
                          suppressMessages(.countBug(relevant.sigs, siglengths)))
  crit_n     <- quantile(sim_vals, 1 - alpha)
  
  # Bootstrap CI on the quantile estimate
  if (!is.null(ci)) {
    boot_quants <- replicate(ci_nsim, {
      quantile(sample(sim_vals, replace = TRUE), 1 - alpha)
    })
    ci_lo <- quantile(boot_quants, (1 - ci) / 2)
    ci_hi <- quantile(boot_quants, 1 - (1 - ci) / 2)
  } else {
    ci_lo <- ci_hi <- NA_real_
  }
  
  list(
    critical_n     = crit_n,
    ci_lower       = ci_lo,
    ci_upper       = ci_hi,
    simulated_values = sim_vals,
    alpha          = alpha,
    quantile_label = paste0(round((1 - alpha) * 100), "th percentile")
  )
}


#' Visualise the Monte Carlo null distribution with threshold and CI
#'
#' Plots a density / histogram of the simulated max-counts from
#' \code{\link{getCriticalN}}, marks the critical threshold, shades the
#' rejection region, and adds a confidence-interval bar for the threshold.
#'
#' @param cn_result   The list returned by \code{\link{getCriticalN}}.
#' @param obs_max     Optional integer.  If supplied, an observed maximum
#'   frequency is marked with a vertical dashed line.
#' @param title       Plot title.
#' @param xlab        X-axis label.
#' @param fill_col    Bar fill colour (default \code{"steelblue"}).
#' @param alpha_fill  Transparency of rejection-region shading (default 0.25).
#'
#' @return A \code{ggplot2} object (invisibly).
#' @export
#'
#' @examples
#' res <- getCriticalN(relevant.sigs, my.siglengths)
#' plotCriticalN(res, obs_max = 7)

plotCriticalN <- function(cn_result,
                          obs_max    = NULL,
                          title      = "Monte Carlo normal distribution of maximum taxon frequency",
                          xlab       = "Maximum taxon count in simulated signatures",
                          fill_col   = "steelblue",
                          alpha_fill = 0.25) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotCriticalN().")
  
  sim_vals <- cn_result$simulated_values
  crit_n   <- cn_result$critical_n
  ci_lo    <- cn_result$ci_lower
  ci_hi    <- cn_result$ci_upper
  alpha    <- cn_result$alpha
  
  df <- data.frame(x = sim_vals)
  
  # Bin width: Sturges-like, at least 1
  bw <- max(1, diff(range(sim_vals)) / (1 + log2(length(sim_vals))))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x)) +
    # Histogram bars
    ggplot2::geom_histogram(
      binwidth = bw,
      fill     = fill_col,
      colour   = "white",
      alpha    = 0.50
    ) +
    # Shaded rejection region (right tail)
    ggplot2::stat_bin(
      binwidth = bw,
      ggplot2::aes(y = ggplot2::after_stat(count),
                   fill = ggplot2::after_stat(x) >= crit_n),
      geom  = "bar",
      colour = "white",
      alpha  = alpha_fill
    ) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "transparent", "TRUE" = "firebrick"),
      guide  = "none"
    ) +
    # Critical-N threshold line
    ggplot2::geom_vline(
      xintercept = crit_n,
      colour     = "firebrick",
      linewidth  = 1.1,
      linetype   = "solid"
    ) +
    ggplot2::annotate(
      "text",
      x      = crit_n,
      y      = Inf,
      label  = paste0("  Critical N = ", round(crit_n),
                      "\n  (", cn_result$quantile_label, ")"),
      hjust  = 0,
      vjust  = 1.4,
      colour = "firebrick",
      size   = 3.5
    )
  
  # Bootstrap CI bar (only if available)
  if (!is.na(ci_lo) && !is.na(ci_hi)) {
    y_ci <- -max(table(cut(sim_vals, breaks = 30))) * 0.06  # just below x-axis
    p <- p +
      ggplot2::annotate(
        "segment",
        x     = ci_lo, xend = ci_hi,
        y     = y_ci,  yend = y_ci,
        colour = "firebrick", linewidth = 1.4
      ) +
      ggplot2::annotate(
        "segment",
        x = ci_lo, xend = ci_lo,
        y = y_ci - abs(y_ci) * 0.4,
        yend = y_ci + abs(y_ci) * 0.4,
        colour = "firebrick", linewidth = 1.0
      ) +
      ggplot2::annotate(
        "segment",
        x = ci_hi, xend = ci_hi,
        y = y_ci - abs(y_ci) * 0.4,
        yend = y_ci + abs(y_ci) * 0.4,
        colour = "firebrick", linewidth = 1.0
      ) +
      ggplot2::annotate(
        "text",
        x = (ci_lo + ci_hi) / 2, y = y_ci,
        label  = paste0("95% CI [", round(ci_lo), ", ", round(ci_hi), "]"),
        vjust  = 2.2,
        colour = "firebrick",
        size   = 3
      )
  }
  
  # Optional observed maximum
  if (!is.null(obs_max)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = obs_max,
        colour     = "darkorange",
        linewidth  = 1,
        linetype   = "dashed"
      ) +
      ggplot2::annotate(
        "text",
        x      = obs_max,
        y      = Inf,
        label  = paste0("  Observed max = ", obs_max),
        hjust  = 0,
        vjust  = 2.8,
        colour = "darkorange",
        size   = 3.5
      )
  }
  
  p <- p +
    ggplot2::labs(
      title    = title,
      subtitle = paste0("n = ", length(sim_vals), " simulations;  ",
                        "shaded region = top ", round(alpha * 100), "% tail"),
      x        = xlab,
      y        = "Count"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40")
    )
  
  print(p)
  invisible(p)
}


significantTaxa <- function(sigs,
                            cn_result,
                            n         = 20,
                            plot      = TRUE,
                            plot_type = c("lollipop", "bar"),
                            title     = "Taxon frequency vs. Monte Carlo threshold") {
  
  plot_type <- match.arg(plot_type)
  
  freq_tab  <- frequencySigs(sigs, n = n)
  crit_n    <- cn_result$critical_n
  sim_vals  <- cn_result$simulated_values
  
  results <- data.frame(
    taxon      = names(freq_tab),
    count      = as.integer(freq_tab),
    critical_n = as.numeric(crit_n),
    stringsAsFactors = FALSE
  )
  results$significant  <- results$count > results$critical_n
  # Empirical p-value: how often does any taxon reach >= observed count by chance?
  results$p_empirical  <- vapply(results$count, function(obs)
    mean(sim_vals >= obs), numeric(1))
  
  # Reorder for plotting (descending frequency, but keep factor order for ggplot)
  results <- results[order(results$count, decreasing = FALSE), ]
  results$taxon <- factor(results$taxon, levels = results$taxon)
  
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
      stop("Package 'ggplot2' is required for plotting in significantTaxa().")
    
    # Colour palette
    sig_colours <- c("TRUE" = "#c0392b", "FALSE" = "#2980b9")
    sig_labels  <- c("TRUE" = paste0("Significant (> Critical N = ",
                                     round(crit_n), ")"),
                     "FALSE" = "Not significant")
    
    if (plot_type == "lollipop") {
      p <- ggplot2::ggplot(results,
                           ggplot2::aes(x     = count,
                                        y     = taxon,
                                        colour = factor(significant))) +
        # Threshold reference line
        ggplot2::geom_vline(
          xintercept = crit_n,
          colour     = "grey30",
          linetype   = "dashed",
          linewidth  = 0.8
        ) +
        # CI shading band around threshold
        {if (!is.na(cn_result$ci_lower))
          ggplot2::annotate("rect",
                            xmin  = cn_result$ci_lower,
                            xmax  = cn_result$ci_upper,
                            ymin  = -Inf, ymax = Inf,
                            fill  = "grey50", alpha = 0.12)
        } +
        # Lollipop stems
        ggplot2::geom_segment(
          ggplot2::aes(x = 0, xend = count, yend = taxon),
          linewidth = 0.7,
          alpha     = 0.6
        ) +
        # Lollipop heads
        ggplot2::geom_point(size = 4) +
        # p-value labels on significant taxa
        ggplot2::geom_text(
          data  = results[results$significant, ],
          ggplot2::aes(label = paste0("p=", signif(p_empirical, 2))),
          hjust = -0.25,
          size  = 3,
          colour = "#c0392b"
        ) +
        ggplot2::scale_colour_manual(
          values = sig_colours,
          labels = sig_labels,
          name   = NULL
        ) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0.12))
        ) +
        ggplot2::annotate(
          "text",
          x     = crit_n,
          y     = nrow(results) + 0.6,
          label = paste0("Critical N = ", round(crit_n)),
          hjust = 0.5,
          vjust = 0,
          size  = 3,
          colour = "grey30"
        )
      
    } else { # bar
      p <- ggplot2::ggplot(results,
                           ggplot2::aes(x    = count,
                                        y    = taxon,
                                        fill = factor(significant))) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_vline(
          xintercept = crit_n,
          colour     = "grey20",
          linetype   = "dashed",
          linewidth  = 0.9
        ) +
        {if (!is.na(cn_result$ci_lower))
          ggplot2::annotate("rect",
                            xmin  = cn_result$ci_lower,
                            xmax  = cn_result$ci_upper,
                            ymin  = -Inf, ymax = Inf,
                            fill  = "grey50", alpha = 0.15)
        } +
        ggplot2::geom_text(
          data  = results[results$significant, ],
          ggplot2::aes(label = paste0("p=", signif(p_empirical, 2))),
          hjust = -0.2,
          size  = 3,
          colour = "white"
        ) +
        ggplot2::scale_fill_manual(
          values = sig_colours,
          labels = sig_labels,
          name   = NULL
        ) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0.10))
        )
    }
    
    p <- p +
      ggplot2::labs(
        title    = title,
        subtitle = paste0("Dashed line = critical N (", cn_result$quantile_label,
                          ");  grey band = 95% CI"),
        x        = "Observed frequency (number of signatures)",
        y        = NULL
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        legend.position  = "bottom",
        plot.title       = ggplot2::element_text(face = "bold"),
        plot.subtitle    = ggplot2::element_text(colour = "grey40"),
        axis.text.y      = ggplot2::element_text(
          face   = ifelse(levels(results$taxon) %in%
                            results$taxon[results$significant],
                          "bold.italic", "plain")),
        panel.grid.major.y = ggplot2::element_blank()
      )
    
    print(p)
  }
  
  invisible(results[order(results$count, decreasing = TRUE), ])
}
