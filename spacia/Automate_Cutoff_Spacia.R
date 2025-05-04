library(diptest)
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

#' Find Exact Correlation Cutoffs for Target Gene Counts
#'
#' This function determines the exact correlation cutoffs that result in 
#' specific numbers of correlated genes. It uses a binary search algorithm 
#' to find precise cutoffs for gene counts ranging from 20 to 2.
#'
#' @param cors Numeric vector of correlation values.
#'
#' @return A data frame with two columns:
#'   \item{cor_cutoffs}{Numeric vector of correlation cutoffs}
#'   \item{num_cor_genes}{Integer vector of corresponding gene counts (20 to 2)}
#'
#' @details
#' The function performs a binary search for each target gene count (20 to 2)
#' to find the correlation cutoff that yields exactly that number of genes.
#' If an exact match is not found, it returns the closest approximation.
#'
#' @examples
#' set.seed(123)
#' cors <- runif(1000, -1, 1)
#' result <- findExactCorCutoffs(cors)
#' print(result)
#'
#' @export
findExactCorCutoffs <- function(cors) {
  # Define target gene counts (from 20 to 2)
  target_gene_counts <- 20:2
  
  # Initialize vector to store correlation cutoffs
  cor_cutoffs <- numeric(length(target_gene_counts))
  
  # Iterate through each target gene count
  for (i in seq_along(target_gene_counts)) {
    target <- target_gene_counts[i]
    
    # Binary search to find the exact cutoffs
    low <- 0
    high <- 1
    
    while (high - low > 1e-6) {  # Precision threshold
      mid <- (low + high) / 2
      count <- sum(abs(cors) > mid)
      
      if (count > target) {
        low <- mid
      } else if (count < target) {
        high <- mid
      } else {
        cor_cutoffs[i] <- mid
        break
      }
    }
    
    # If we didn't find an exact match, use the closest approximation
    if (cor_cutoffs[i] == 0) {
      cor_cutoffs[i] <- (low + high) / 2
    }
  }
  
  # Create a data frame with results
  cor_df <- data.frame(
    cor_cutoffs = cor_cutoffs,
    num_cor_genes = target_gene_counts
  )
  
  return(cor_df)
}




#' Check for Bimodality and Find Midpoint
#'
#' This function analyzes a numeric vector to determine if its distribution
#' is bimodal and, if so, calculates the midpoint between the two modes.
#'
#' @param signature A numeric vector representing the distribution to be analyzed.
#' @param alpha Numeric. Significance level for Hartigan's dip test (default: 0.05).
#' @param min_peak_ratio Numeric. Minimum ratio of the lower peak to the higher peak (default: 0.3).
#' @param min_valley_ratio Numeric. Minimum ratio of valley depth to lower peak height (default: 0.5).
#'
#' @return A list with the following components:
#'   \item{is_bimodal}{Logical. TRUE if the distribution is determined to be bimodal, FALSE otherwise.}
#'   \item{midpoint}{Numeric. The x-value of the deepest valley between the two highest peaks if bimodal, NULL otherwise.}
#'   \item{dip_pvalue}{Numeric. P-value from Hartigan's dip test.}
#'   \item{peak_ratio}{Numeric. Ratio of the lower peak to the higher peak.}
#'   \item{valley_ratio}{Numeric. Ratio of valley depth to lower peak height.}
#'
#' @examples
#' set.seed(123)
#' bimodal_data <- c(rnorm(500, mean = -2), rnorm(500, mean = 2))
#' result <- check_bimodality_and_midpoint(bimodal_data)
#' print(result)
#'
#' @import diptest
#' @importFrom stats density
#'
#' @export
check_bimodality_and_midpoint <- function(signature, alpha = 0.05, min_peak_ratio = 0.05, min_valley_ratio = 0.1) {
  # Perform Hartigan's dip test, suppressing messages
  dip_result <- suppressMessages(diptest::dip.test(signature))
  
  # If dip test doesn't suggest multimodality, return early
  if (dip_result$p.value >= alpha) {
    return(list(is_bimodal = FALSE, midpoint = NULL))
  }
  
  # Estimate density
  density_result <- stats::density(signature)
  
  # Find peaks and valleys
  peaks <- findPeaks(density_result$y)
  valleys <- findValleys(density_result$y)
  
  if (length(peaks) >= 2) {
    # Get the two highest peaks
    top_two_peaks <- peaks[order(density_result$y[peaks], decreasing = TRUE)[1:2]]
    peak_heights <- density_result$y[top_two_peaks]
    peak_ratio <- min(peak_heights) / max(peak_heights)
    
    if (peak_ratio > min_peak_ratio) {
      # Find the deepest valley between the two highest peaks
      valley_between <- valleys[valleys > min(top_two_peaks) & valleys < max(top_two_peaks)]
      
      if (length(valley_between) > 0) {
        deepest_valley <- valley_between[which.min(density_result$y[valley_between])]
        valley_depth <- min(peak_heights) - density_result$y[deepest_valley]
        valley_ratio <- valley_depth / min(peak_heights)
        
        if (valley_ratio > min_valley_ratio) {
          midpoint <- density_result$x[deepest_valley]
          return(list(is_bimodal = TRUE, 
                      midpoint = midpoint, 
                      dip_pvalue = dip_result$p.value, 
                      peak_ratio = peak_ratio, 
                      valley_ratio = valley_ratio))
        }
      }
    }
  }
  
  return(list(is_bimodal = FALSE, midpoint = NULL))
}




#' Find Peaks in a Numeric Vector
#'
#' @param x Numeric vector
#' @return Integer vector of peak indices
#' @keywords internal
findPeaks <- function(x) which(diff(sign(diff(x))) < 0) + 1

#' Find Valleys in a Numeric Vector
#'
#' @param x Numeric vector
#' @return Integer vector of valley indices
#' @keywords internal
findValleys <- function(x) which(diff(sign(diff(x))) > 0) + 1




#' Plot Density Distribution with Midpoint
#'
#' This function creates a density plot of the signature distribution,
#' highlighting the midpoint and providing relevant statistical information.
#' It also saves the plot as a PDF file in a specified folder.
#'
#' @param signature Numeric vector. The signature values to plot.
#' @param midpoint Numeric. The calculated midpoint to highlight on the plot.
#' @param cor_cutoff Numeric. The correlation cutoff used in the analysis.
#' @param quantile Numeric. The calculated quantile value.
#' @param num_cor_genes Integer. The number of correlated genes used.
#' @param is_bimodal Logical. Whether the distribution is bimodal.
#' @param receiving_gene Character. The name of the receiving gene.
#' @param save_folder Character. The folder path to save the PDF plot.
#'
#' @return A ggplot object of the density plot.
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @examples
#' # Assuming necessary variables are defined
#' plot <- plot_density_with_midpoint(signature, midpoint, cor_cutoff, quantile, 
#'                                    num_cor_genes, is_bimodal, "GENE1", "plots/")
#' print(plot)
#'
#' @export
plot_density_with_midpoint <- function(signature, midpoint, cor_cutoff, quantile, 
                                       num_cor_genes, is_bimodal, receiving_gene, 
                                       save_folder) {
  # Construct the plot title
  plot_title <- sprintf("Signature Distribution of %s (Cutoff: %.4f, Quantile: %.4f, Genes: %d)",
                        receiving_gene, cor_cutoff, quantile, num_cor_genes)
  
  # Modify title if distribution is not bimodal
  if (!is_bimodal) {
    plot_title <- paste(plot_title, "(Non-bimodal, using median)")
  }
  
  # Create the ggplot object
  plot <- ggplot(data.frame(signature = signature), aes(x = signature)) +
    geom_density(fill = "lightblue", alpha = 0.7) +
    geom_vline(xintercept = midpoint, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = plot_title,
         x = "Signature Value",
         y = "Density") +
    theme_minimal() +
    annotate("text", x = midpoint, y = 0, label = sprintf("Midpoint: %.2f", midpoint),
             vjust = -0.5, hjust = -0.1, color = "red")
  
  # Save the plot as a PDF file
  pdf_file <- file.path(save_folder, paste0(receiving_gene,"_automated_cutoff", ".pdf"))
  pdf(pdf_file, width = 10, height = 7)  # Adjust width and height as needed
  print(plot)
  dev.off()
  
  # Return the plot object
  return(plot)
}




#' Automated Cutoff Generator for Spacia's bag binarization
#'
#' This function automates the process of finding optimal correlation cutoffs
#' and quantile cutoffs for Spacia R's receiving_gene binarization. 
#' It performs correlation analysis, searches for bimodal distributions, 
#' and generates visualizations of the results.
#'
#' @param receiving_gene Character. The name of the receiving gene to analyze.
#' @param counts_receiver data.table. Contains gene expression counts.
#'        The first column should be sample IDs, and other columns should be gene counts.
#' @param save_folder Character. The path to the folder where plots will be saved.
#' @param exp_sender List. Bags created for each receiving cell (used in signature calculation).
#' @param backup_num_cor_genes Integer. Number of correlated genes to use if no bimodal 
#'        distribution is found. Should be between 2 and 20. Default is 10.
#' @param alpha Numeric. Significance level for Hartigan's dip test. Default is 0.05.
#' @param min_peak_ratio Numeric. Minimum ratio of the lower peak to the higher peak. Default is 0.05.
#' @param min_valley_ratio Numeric. Minimum ratio of valley depth to lower peak height. Default is 0.1.
#'
#' @return A list containing the results of the analysis, including:
#'   \item{data}{A data frame with the following columns:
#'     \itemize{
#'       \item cor_cutoff: The optimal correlation cutoff found
#'       \item midpoint: The midpoint of the bimodal distribution or median
#'       \item quantile: The quantile corresponding to the midpoint
#'       \item num_cor_genes: The number of correlated genes at the cutoff
#'       \item is_bimodal: Logical indicating if a bimodal distribution was found
#'     }
#'   }
#'   \item{plot}{The ggplot object of the density plot}
#'
#' @import data.table
#' @importFrom stats cor ecdf median
#' @importFrom utils head
#'
#' @examples
#' # Assuming necessary data and functions are loaded
#' result <- AutomatedCutoffGenerator("GENE1", counts_data, "plots/", exp_sender_data)
#' print(result$data)
#' print(result$plot)
#'
#' @export
AutomatedCutoffGenerator <- function(receiving_gene, counts_receiver, save_folder, exp_sender, 
                                     backup_num_cor_genes = 10, alpha = 0.05, 
                                     min_peak_ratio = 0.05, min_valley_ratio = 0.1) {
  
  # Input validation
  if (!is.character(receiving_gene) || length(receiving_gene) != 1) {
    stop("receiving_gene must be a single character string")
  }
  if (!data.table::is.data.table(counts_receiver)) {
    stop("counts_receiver must be a data.table")
  }
  if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
  }
  if (!receiving_gene %in% colnames(counts_receiver)) {
    stop("receiving_gene not found in counts_receiver")
  }
  if (!is.numeric(backup_num_cor_genes) || backup_num_cor_genes < 2 || backup_num_cor_genes > 20) {
    stop("backup_num_cor_genes must be an integer between 2 and 20")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }
  if (!is.numeric(min_peak_ratio) || min_peak_ratio <= 0 || min_peak_ratio >= 1) {
    stop("min_peak_ratio must be a numeric value between 0 and 1")
  }
  if (!is.numeric(min_valley_ratio) || min_valley_ratio <= 0 || min_valley_ratio >= 1) {
    stop("min_valley_ratio must be a numeric value between 0 and 1")
  }
  
  # Calculate correlations
  cors <- cor(counts_receiver[[receiving_gene]], 
              counts_receiver[, -1, with = FALSE], use = "pairwise.complete.obs")[1,]
  
  # Find exact correlation cutoffs
  cor_df <- findExactCorCutoffs(cors)
  
  result <- NULL
  backup_iteration <- NULL
  plot <- NULL
  
  for (i in 1:nrow(cor_df)) {
    cor_cutoff <- cor_df$cor_cutoffs[i]
    num_cor_genes <- cor_df$num_cor_genes[i]
    
    keep <- abs(cors) > cor_cutoff
    selected_genes <- names(cors)[keep]
    
    # Signature calculation
    signature <- counts_receiver[
      get(names(counts_receiver)[1]) %in% names(exp_sender),
      .SD,
      .SDcols = selected_genes
    ][, {
      result <- colMeans(t(as.matrix(.SD)) * cors[keep])
      .(signature = result)
    }]
    
    # Check bimodality and find midpoint with new parameters
    bimodal_check <- check_bimodality_and_midpoint(signature$signature, alpha, min_peak_ratio, min_valley_ratio)
    
    if (bimodal_check$is_bimodal) {
      midpoint <- bimodal_check$midpoint
      quantile_result <- ecdf(signature$signature)(midpoint)
      
      result <- data.frame(
        cor_cutoff = cor_cutoff,
        midpoint = midpoint,
        quantile = quantile_result,
        num_cor_genes = num_cor_genes,
        is_bimodal = TRUE
      )
      
      # Generate the density plot
      plot <- plot_density_with_midpoint(signature$signature, 
                                         midpoint, 
                                         cor_cutoff, 
                                         quantile_result, 
                                         num_cor_genes, 
                                         TRUE,
                                         receiving_gene,
                                         save_folder)
      
      break  # Exit the loop once we find a bimodal distribution
    }
    
    # Store the backup iteration data for potential use later
    if (num_cor_genes == backup_num_cor_genes) {
      backup_iteration <- list(
        signature = signature$signature,
        cor_cutoff = cor_cutoff,
        num_cor_genes = num_cor_genes
      )
    }
  }
  
  # If no bimodal distribution was found, use the backup cutoff and median
  if (is.null(result)) {
    if (is.null(backup_iteration)) {
      stop("No valid backup iteration found. Check your backup_num_cor_genes value.")
    }
    
    signature <- backup_iteration$signature
    cor_cutoff <- backup_iteration$cor_cutoff
    num_cor_genes <- backup_iteration$num_cor_genes
    
    midpoint <- median(signature)
    quantile_result <- 0.5  # by definition, the median is at the 0.5 quantile
    
    result <- data.frame(
      cor_cutoff = cor_cutoff,
      midpoint = midpoint,
      quantile = quantile_result,
      num_cor_genes = num_cor_genes,
      is_bimodal = FALSE
    )
    
    # Generate the density plot
    plot <- plot_density_with_midpoint(signature, 
                                       midpoint, 
                                       cor_cutoff, 
                                       quantile_result, 
                                       num_cor_genes, 
                                       FALSE,
                                       receiving_gene,
                                       save_folder)
    
    cat(sprintf("No bimodal distribution found. Using backup cutoff for %s (num_cor_genes = %d) and median.\n", 
                receiving_gene, 
                backup_num_cor_genes))
  }
  
  # Return both the result data frame and the plot
  return(list(data = result, plot = plot))
}



#' Generate Plots for Multiple Correlation Cutoffs with Visible Quantile Labels
#'
#' This function generates density plots for gene signatures at different correlation cutoffs,
#' allowing users to visualize the distribution and manually annotate suitable cutoffs.
#' It includes clearly visible labeled quantile lines for easy reference.
#'
#' @param receiving_gene Character. The name of the receiving gene to analyze.
#' @param counts_receiver data.table. Contains gene expression counts.
#'        The first column should be sample IDs, and other columns should be gene counts.
#' @param save_folder Character. The path to the folder where plots will be saved.
#' @param exp_sender List. Bags created for each receiving cell (used in signature calculation).
#'
#' @return A grob object containing the combined plot of all cutoffs.
#'
#' @import ggplot2
#' @import gridExtra
#' @import RColorBrewer
#' @importFrom stats cor quantile
#' @importFrom grDevices pdf dev.off
CutoffPlotGenerator <- function(receiving_gene, counts_receiver, output_file, exp_sender) {
  # Input validation
  if (!is.character(receiving_gene) || length(receiving_gene) != 1) {
    stop("receiving_gene must be a single character string")
  }
  if (!data.table::is.data.table(counts_receiver)) {
    stop("counts_receiver must be a data.table")
  }
  if (!receiving_gene %in% colnames(counts_receiver)) {
    stop("receiving_gene not found in counts_receiver")
  }
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Calculate correlations
  cors <- cor(counts_receiver[[receiving_gene]], 
              counts_receiver[, -1, with = FALSE], use = "pairwise.complete.obs")[1,]
  
  # Find exact correlation cutoffs
  cor_df <- findExactCorCutoffs(cors)
  
  # Define quantiles and colors
  quantiles <- seq(0.1, 0.9, by = 0.1)
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
              "#FFFF33", "#A65628", "#F781BF", "#999999")
  
  # Create a list to store all plots
  all_plots <- list()
  
  # Generate plots for each correlation cutoff
  for (i in 1:nrow(cor_df)) {
    cor_cutoff <- cor_df$cor_cutoffs[i]
    num_cor_genes <- cor_df$num_cor_genes[i]
    
    # Select genes based on correlation cutoff
    keep <- abs(cors) > cor_cutoff
    selected_genes <- names(cors)[keep]
    
    # Calculate signature
    signature <- counts_receiver[
      get(names(counts_receiver)[1]) %in% names(exp_sender),
      .SD,
      .SDcols = selected_genes
    ][, {
      result <- colMeans(t(as.matrix(.SD)) * cors[keep])
      .(signature = result)
    }]
    
    # Calculate quantile values
    quantile_values <- quantile(signature$signature, probs = quantiles)
    
    # Create the density plot
    plot <- ggplot(data.frame(signature = signature$signature), aes(x = signature)) +
      geom_density(fill = "lightblue", alpha = 0.7) +
      labs(title = sprintf("%s (Cutoff: %.4f, Genes: %d)",
                           receiving_gene, cor_cutoff, num_cor_genes),
           x = "Signature Value",
           y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10),
            plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt"))
    
    # Get the y-range of the plot for label positioning
    y_range <- ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range
    
    # Add quantile lines and labels to the plot
    for (j in 1:length(quantiles)) {
      plot <- plot +
        geom_vline(xintercept = quantile_values[j], 
                   color = colors[j], 
                   linetype = "dashed", 
                   linewidth = 0.5) +
        annotate("text", 
                 x = quantile_values[j], 
                 y = y_range[2] * 1.05,
                 label = sprintf("%.1f", quantiles[j]), 
                 color = colors[j], 
                 angle = 90, 
                 size = 2.5,
                 vjust = 0)
    }
    
    all_plots[[i]] <- plot
  }
  
  # Create a legend plot
  legend_data <- data.frame(
    x = 1:9,
    y = rep(1, 9),
    label = sprintf("%.1f", quantiles)
  )
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = factor(x))) +
    geom_point() +
    scale_color_manual(values = colors, labels = legend_data$label) +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(color = "Quantile") +
    guides(color = guide_legend(nrow = 3, byrow = TRUE))
  
  # Arrange plots in a grid
  n_plots <- length(all_plots)
  n_cols <- 5
  n_rows <- ceiling(n_plots / n_cols)
  
  grid_plots <- c(all_plots, list(legend_plot))
  if (length(grid_plots) < n_cols * n_rows) {
    grid_plots <- c(grid_plots, replicate(n_cols * n_rows - length(grid_plots), ggplot() + theme_void()))
  }
  
  # Combine plots
  combined_plot <- gridExtra::arrangeGrob(grobs = grid_plots, ncol = n_cols)
  
  # Save the combined plot as a PDF file
  ggsave(output_file, combined_plot, width = 20, height = 4 * n_rows, limitsize = FALSE)
}