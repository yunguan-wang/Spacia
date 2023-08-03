suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(ggmcmc))

#' Create density, trace, autocorrelation, and PSRF plots 
#' for b and betas and return the ggmcmc::ggs() object, S
#'
#' @param beta_matrix matrix of the beta, where column names are populated
#' @param b_matrix matrix of the b, where column names are populated
#' @param nwarm number of warm up iterations (burn-in periods)
#' @param ntotal total number of MCMC iterations
#' @param nthin number used to thin and store posteriors
#' @param nchain number of chains used in MCMC
#' @param job_id string that specifies the job ID
#' @param output_path string of output path
#' @param ext extension used for the ggsave. Device to use. 
#'   Can either be a device function (e.g. png), or one of 
#'   "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", 
#'   "bmp", "svg" or "wmf" (windows only).
#'
#' @return ggmcmc::ggs() object
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' S <- BetaB2MCMCPlots(
#'     beta_matrix, 
#'     b_matrix, 
#'     20000, 
#'     50000, 
#'     100, 
#'     2,
#'     "test",
#'     "./output_dir/",
#'     "pdf")
#'     
#' }
BetaB2MCMCPlots <- function(beta_matrix, 
                            b_matrix, 
                            nwarm, 
                            ntotal, 
                            nthin, 
                            nchain,
                            job_id,
                            output_path,
                            ext 
                            ){
  
  # Create variables for the coda::mcmc object
  start = ntotal - nwarm + nthin
  end = ntotal
  niter = ntotal - nwarm
  
  # number of rows saved for each chain
  nsaved = 1 + (1 + floor((niter - 1) / nthin))
  
  # for mcmc.list() to take in
  list_of_mcmcs = vector(mode = "list", length = nchain)
  
  # Populate list_of_mcmcs by cbind b and beta
  for (chain in 1:nchain) {
    start_idx = ((chain-1)*nsaved+1) + 1 # We don't want iter 1
    end_idx = chain*nsaved
    
    data = cbind(
      beta_matrix[start_idx:end_idx,],
      b_matrix[start_idx:end_idx,2,drop=F]
    )
    
    list_of_mcmcs[[chain]] = mcmc(
      data = data,
      start = start,
      end = end,
      thin = nthin
    )
  }
  
  # Make mcmc.list obj
  mcmc_list_obj <- mcmc.list(list_of_mcmcs)
  
  # coda::mcmc.list to ggmcmc::ggs
  S = ggs(mcmc_list_obj)
  
  # Reorder levels
  nbeta = dim(beta_matrix)[2]
  S$Parameter <- factor(S$Parameter, levels = c(c("b.2"), paste("beta.", 1:nbeta, sep="")))
  
  # Calculate the appropriate width and height of the plot
  # When ncol(beta) = 50, (width, height) = (30, 30), (30, 30), (30, 90) cm works 
  size = c(0, 0, 0)
  if (nbeta <= 4){ # one row
    size[1] = 3
    size[2] = round(nbeta/50 * 90)
  } else if (nbeta >= 66) { # max rows, ggsave supports max 127cm
    size[1] = 120
    size[2] = round(nbeta/50 * 90)
  } else { # calculate appropriate rows
    size[1] = round(nbeta/50 * 30)
    size[2] = round(nbeta/50 * 90)
  }
  
  size[3] = nchain * 15
  
  # Plot density, trace plot, and autocorrelation
  ggs_density(S)+ facet_wrap(~ Parameter, ncol = 5, scales = "free")
  ggsave(paste(output_path, job_id,'_density.', ext, sep=''), 
         width = 30, height = size[1], units = "cm")
  
  ggs_traceplot(S) + facet_wrap(~ Parameter, ncol = 5, scales = "free") 
  ggsave(paste(output_path, job_id,'_trace.', ext, sep=''), 
         width = 30, height = size[1], units = "cm")
  
  ggs_autocorrelation(S)
  ggsave(paste(output_path, job_id,'_autocorr.', ext, sep=''), 
         width = size[3], height = size[2], units = "cm")
  
  
  if ( nchain <= 1 ) {
    # We need multiple chains to assess PSRF
    warning("The number of chains is less than 2. Skipping PSRF plot.")
  } else {
    ggs_Rhat(S) + xlab("R_hat")
    ggsave(paste(output_path, job_id,'_psrf.', ext, sep=''), 
           width = round(size[3]/3), height = round(size[2]/6), units = "cm")
  }

  
  # return the ggs object since we may want more plots
  return(S)
  
}
