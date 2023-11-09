# library(coda)

#### Get hyperparameters ####
getHyperPars <- function(tidydata){
  D = tidydata$nfeature_inst # pos (1) + exp (D-1)
  res = list(
    hp_mu_beta=rep(0, D),
    hp_mu_b=rep(0,2),
    
    hp_Sig_beta = diag(c(1,rep(1,D-1)),D), 
    hp_Sig_b = diag(c(1,rep(1,1)),2)

    # Testing model priors

    # hp_Sig_beta = diag(c(5, rep(5, D-1)),D), 
    # hp_Sig_b = diag(c(5, rep(5, 1)),2)

    # hp_Sig_beta = diag(c(1000,rep(1000,D-1)),D), 
    # hp_Sig_b = diag(c(1000,rep(1000,1)),2)

    # hp_Sig_beta = diag(c(100,rep(100,D-1)),D), 
    # hp_Sig_b = diag(c(100,rep(100,1)),2)

    # hp_Sig_beta = diag(c(10,rep(10,D-1)),D), 
    # hp_Sig_b = diag(c(10,rep(10,1)),2)
  )
  return(res)
}

#### Get initial values ####
getInits <- function(tidydata,hyperpars){
  beta = unlist(lapply(hyperpars$hp_mu_beta, function(mu_beta) rnorm(1, mu_beta, 10)))
  b = unlist(lapply(hyperpars$hp_mu_b, function(mu_b) rnorm(1, mu_b, 10)))
  delta = unlist(lapply(1:tidydata$nsample,function(i){rbinom(tidydata$ninst[i],1,mean(tidydata$label))}))
  
  res = list(
    beta = beta,
    b = b,
    delta = delta
  )
  return(res)
}


#### Summarize all input data and parameters for mcmc chain #### 
getInputPars <- function(tidydata){
  
  list_hyperpars <- getHyperPars(tidydata)
  list_inits <- getInits(tidydata,list_hyperpars)
  
  tmp=Reduce(rbind, tidydata$feature_inst)
  res = list(
    ## data
    n = tidydata$nsample, # number of bags
    d = tidydata$nfeature_inst, # number of features, pos (1) + exp (d-1)
    m = tidydata$ninst, # number of instances per bag
    membership = tidydata$membership, # membership for instances
    y = tidydata$label, # bag labels
    X1 = cbind(rep(1,dim(tmp)[1]),tmp), # design matrix
    ## hyperparameters 
    hp_mu_beta = list_hyperpars$hp_mu_beta,
    hp_mu_b = list_hyperpars$hp_mu_b,
    hp_Sig_beta = list_hyperpars$hp_Sig_beta,
    hp_Sig_b = list_hyperpars$hp_Sig_b,
    ## model parameters
    beta = list_inits$beta,
    b = list_inits$b,
    delta = list_inits$delta
  )
  return(res)
}


#### Fitting BMIR2 model ####

#### 1 Gibbs iteration in Rcpp ####
MICProB_sampler<-function(tidytrain,
                        tidytest,
                        ntotal,
                        nwarm,
                        nthin,
                        nchain,
                        #scale,
                        return_delta,
                        prior = 1){
  
  cat("=============================================================\n")
  cat(sprintf("Probit Bayesian Multiple Instance Classification\n"))
  
  res_mcmc <- vector("list", nchain)
  
  for(nc in 1:nchain){
    
    # begin time
    start_time <- Sys.time()
    
    parlist <- getInputPars(tidytrain)
    
    y<-parlist$y
    n<-parlist$n
    d<-parlist$d
    m<-parlist$m
    N<-sum(m)
    membership<-parlist$membership
    
    hp_mu_beta<-parlist$hp_mu_beta
    hp_mu_b<-parlist$hp_mu_b
    hp_Sig_beta<-parlist$hp_Sig_beta
    hp_Sig_b<-parlist$hp_Sig_b

    if (prior != 1) {
      hp_Sig_beta = diag(c(prior, rep(prior, d-1)),d)
      hp_Sig_b = diag(c(prior, rep(prior, 1)),2)
      cat(sprintf("prior b and beta resetted.\n"))
    }
    
    beta<-parlist$beta
    b<-parlist$b
    delta<-parlist$delta
    u = rep(0,N)
    z = rep(0,n)
    
    hp_Sig_beta_inv<-solve(hp_Sig_beta)
    hp_Sig_b_inv<-solve(hp_Sig_b)
    
    # posterior variance of b
    X1 <- parlist$X1
    V_b <- solve(hp_Sig_b_inv + crossprod(X1[,1:2], X1[,1:2]))
    
    #cat("=============================================================\n")
    #cat(sprintf("Bayesian Multiple Instance Regression: chain" ,nc, " \n"))
    
    # begin time
    # start_time <- Sys.time()
    
    tick = 0.2
    
    # Gibbs sampling (warming up)
    
    cat("=============================================================\n")
    cat("Start warming up",nwarm,"MCMC samples!\n")
    cat("Progress: ")
    
    for(iter in 1:nwarm){
      if(iter %in% seq(round(tick*nwarm),nwarm,by=round(tick*nwarm))){
        cat(100*iter/nwarm,"% ...")
      }
      mcmc_res <- MICProB_1Gibbs_cpp(Xb = X1[,2,drop=F],Xbeta=X1[,-c(1,2), drop = F],
                                     y = y,
                                      ninst = m,
                                      hp_mu_beta = hp_mu_beta,
                                      hp_mu_b,
                                      hp_Sig_beta,
                                      hp_Sig_b,
                                      beta,
                                      b,
                                      delta,
                                      u,
                                      z,
                                      hp_Sig_beta_inv,
                                      hp_Sig_b_inv,
                                      V_b)

      # update parameters
      beta = mcmc_res$beta
      b = mcmc_res$b
      delta = mcmc_res$delta
      u = mcmc_res$u
      z = mcmc_res$z
      
    } # end warm-up
    cat("\n")
    cat("Finish warming up!\n")
    cat("-------------------------------------------------------------\n")
    
    niter = ntotal - nwarm
    nsave = 1 + floor((niter - 1) /nthin)
    
    # posterior quantities to be saved
    beta_post<-matrix(NA,nrow=nsave,ncol=length(beta))
    b_post<-matrix(NA,nrow=nsave,ncol=length(b))
    delta_post<-matrix(NA,nrow=nsave,ncol=length(delta))
    
    pip_1chain<-rep(0,length(delta))
    mcmc_1chain <- list()
    
    cat("Start extracting",niter,"MCMC samples!\n")
    cat("Progress :")
    for(iter in 1:niter){
      if(iter %in% seq(round(tick*niter),niter,by=round(tick*niter))){
        cat(100*iter/niter,"% ...")
      }
      mcmc_res <- MICProB_1Gibbs_cpp(Xb = X1[,2,drop=F],Xbeta=X1[,-c(1,2), drop = F],
                                     y = y,
                                      ninst = m,
                                      hp_mu_beta = hp_mu_beta,
                                      hp_mu_b,
                                      hp_Sig_beta,
                                      hp_Sig_b,
                                      beta,
                                      b,
                                      delta,
                                      u,
                                      z,
                                      hp_Sig_beta_inv,
                                      hp_Sig_b_inv,
                                      V_b)

      # update parmaeters
      beta = mcmc_res$beta
      b = mcmc_res$b
      delta = mcmc_res$delta
      u = mcmc_res$u
      z = mcmc_res$z
      
      # save posterior samples
      if(iter %in% seq(nthin,niter,by=nthin)){ # thinning delta
        if(return_delta){
          delta_post[iter/nthin,]<-delta
        }
        pip_1chain = pip_1chain + delta
        beta_post[iter/nthin,]<-beta
        b_post[iter/nthin,]<-b
      }
      
    } # end extracting posterior samples
    
    pip_1chain = pip_1chain / nsave
    cat("\n")
    cat("Finish MCMC sampling!\n")
    cat("=============================================================\n")
    
    # elapsed time
    cat(sprintf("Elapsed time for chain%d=%.3f mins: MCMC sampling is done!\n", nc, difftime(Sys.time(), start_time, units = "mins")))
    
    # output
    mcmc_1chain[["beta"]]<-rbind(parlist$beta,beta_post)
    mcmc_1chain[["b"]]<-rbind(parlist$b,b_post)
    mcmc_1chain[["pip"]]<-pip_1chain
    
    if(return_delta){
      mcmc_1chain[["delta"]]<-rbind(parlist$delta, delta_post)
    } else{
      mcmc_1chain[["delta"]]<-NULL
    }

    res_mcmc[[nc]]<-mcmc_1chain
    
  }
  
  return(res_mcmc)
}
