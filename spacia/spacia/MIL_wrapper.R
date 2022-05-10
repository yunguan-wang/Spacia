########  simulate some very simple data   ###################
# for the purpose of debug and show the format of the input data

# # bag level, expression of one gene in the receiver cell
# # dichotomized by the users
# exp_receiver=c(1,1,0,0)
# 
# # instance features, distances of the sending cells to receiver cells
# # each bag should have at least one instance
# pos_sender=list(
#  c(1,4,3),
#  c(1),
#  c(1,5,2,3),
#  c(3)
# )
# 
# # instance features, expression of genes/pathways in the sending cells
# # rows are genes/pathways, cols are instances
# exp_sender=list(
#  matrix(runif(3*3),ncol=3),
#  matrix(runif(1*3),ncol=3),
#  matrix(runif(4*3),ncol=3),
#  matrix(runif(1*3),ncol=3)
# )
# 
# # MCMC parameters
# ntotal=100
# nwarm=10
# nthin=10
# nchain=1
# 
# # cutoffs on the pi_hats
# thetas=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
# res = MIL_C2Cinter(
#     exp_receiver, pos_sender, exp_sender, ntotal, nwarm, nthin, nchain, thetas)

########3  MIL wrapper  #####################

MIL_C2Cinter<-function(exp_receiver,pos_sender,exp_sender,
  ntotal,nwarm,nthin,nchain,thetas)
{
  # organize into Danyi's original format
  tidy_train=list()
  tidy_train$membership=unlist(sapply(1:length(exp_sender),
    function(i) rep(i,dim(exp_sender[[i]])[1])))           
  tidy_train$ninst=sapply(exp_sender,function(x) dim(x)[1])
  
  tidy_train$label=exp_receiver
  tidy_train$nsample=length(exp_receiver)

  exp_pos_sender=exp_sender
  for (i in 1:length(exp_sender))
    {exp_pos_sender[[i]]=cbind(pos_sender[[i]],exp_pos_sender[[i]])}
  
  tidy_train$feature_inst=exp_pos_sender
  tidy_train$nfeature_inst=dim(exp_pos_sender[[1]])[2]
  
  tidytrain=tidy_train
  tidydata=tidy_train
  
  # invoke Danyi's functions
  res_mcmc=MICProB_sampler(tidytrain,
                            tidytest=NULL,
                            ntotal,
                            nwarm,
                            nthin,
                            nchain,
                            #scale,
                            return_delta=TRUE)
  
  # organize results
  pip=c() # col=nchain, row=number of senders
  for (nc in 1:nchain) {pip=cbind(pip,res_mcmc[[nc]]$pip)}

  b=c() # all samples of all MCMC chains
  for (nc in 1:nchain) {b=rbind(b,res_mcmc[[nc]]$b)}
  
  beta=c() # col=number of features, row=all samples of all chains
  for (nc in 1:nchain) {beta=rbind(beta,res_mcmc[[nc]]$beta[,-1,drop=F])}
  
  FDRs=sapply(thetas,function(theta) { # B-FDRs
    sum((pip>theta)*(1-pip))/sum(pip>theta)
  })
  
  # recalculate a point-estiate of the probability of primary instances 
  # to avoid the negative impact of randomness in MCMC sampling
  # this still is not ideal
  # good enough to use it as a "score" to indicate primary instance
  # but cannot interpret it as a probability (as it should be)
  
  pip_recal=c()
  b0=mean(b[,1])
  b1=mean(b[,2])
  for (i in 1:length(pos_sender))
    {pip_recal=c(pip_recal,pnorm(b0+pos_sender[[i]]*b1))}
  
  # (1) pip is the probabilities of the sender cells being truly responsible for the
  # receiving cells' phenotypes
  # (2) b is the coef for explaining how the distances between sender and receiver cells
  # affect the probabilities of the sender cells being responsible
  # (3) beta are the coefs for explaining the effect of the expression of the genes 
  # in the sender cells
  # (4) FDRs are the Bayes FDRs calculated according to the theta
  # cutoffs given by the users (for defining primary instances)
  # (5) recalculated pip
  return(list(pip=pip,b=b,beta=beta,FDRs=FDRs,pip_recal=pip_recal))
}
