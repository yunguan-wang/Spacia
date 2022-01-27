##########  set Rstudio environment  ############

# this is a problem only specific to the Rstudio of biohpc
Sys.setenv(LIBRARY_PATH = "/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mkl/lib/intel64:/cm/shared/apps/java/oracle/jdk1.7.0_51/lib:/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/compiler/lib/intel64:/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mpi/intel64/lib:/cm/shared/apps/openmpi/gcc/64/2.1.5/lib64:/cm/shared/apps/gcc/5.4.0/lib:/cm/shared/apps/gcc/5.4.0/lib64:/cm/shared/apps/slurm/16.05.8/lib64/slurm:/cm/shared/apps/slurm/16.05.8/lib64")

########  tuning parameters of the simulation  ###########

# number of receiver cells
numbers_receiver=1000

# maximum sender cells for one receiver cell
max_sender=20

# position cutoff to define primary sender cells
# between 0.1 and 2, smaller values means the 
# sender and receiver cells need to be closer to be interacting
cutoff=1

# number of features in the expressio matrix of sender cells
# features could be single genes, but are more likely 
# aggregates of genes (pathways)
nfeature=50

# quantile for dichotomizing the expression of 
# the gene of interest in the receiver cells into binary status
q_receiver=0.5

###########  generate simulation data  ###########

# number of sender cells
numbers_sender=round(runif(numbers_receiver)*(max_sender-1)+1)
# number of primary sender cells
numbers_primary_sender=sapply(numbers_sender,
  function(i) round(runif(1,1,i)))

# function to simulate position
# primary sender cells are closer than non-primary ones
simulate_pos<-function(primary,n,cutoff)
{
  if (n==0) {return(c())}
  if (primary) {runif(n,min=0.1,max=cutoff)} else  {runif(n,min=cutoff,max=2)}
}

# simulate position
pos_sender=sapply(1:length(numbers_sender), function(i) {
  c(simulate_pos(T,numbers_primary_sender[i],cutoff),
    simulate_pos(F,numbers_sender[i]-numbers_primary_sender[i],cutoff))
})

# gold-standard for primary instance status
pip=unlist(pos_sender)<cutoff

# expression data of sender cells
# important: need to keep the expression data non-negative here
# to mimic x or log(x+1), x=raw/normalized counts
exp_sender=sapply(1:numbers_receiver, function(i) {
  matrix(runif(numbers_sender[i]*nfeature),ncol=nfeature)
})

# gold standard coef 
# for regression expression of senders onto the status of receivers
beta=matrix(rnorm(nfeature),ncol=1)

# functional status of the receiver cells
# note: our model investigates all genes/pathways in the sending 
# cells at the same time, but can only investigate
# one gene/pathway in the receiver cells one at a time
exp_receiver=unlist(sapply(1:numbers_receiver, function(i) {
  tmp=exp_sender[[i]] %*% beta
  sum(tmp[1:numbers_primary_sender[i],,drop=F])
}))
exp_receiver=1*(exp_receiver>quantile(exp_receiver,q_receiver))

################  running the MIL model  ###################

# other input parameters
ntotal=8000
nwarm=4000
nthin=10
nchain=1

thetas=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# run MIL
source('~/projects/MIL4ST/MIL_wrapper.R')

Sys.time()
results=MIL_C2Cinter(exp_receiver,pos_sender,exp_sender,
  ntotal,nwarm,nthin,nchain,thetas)
Sys.time()

###########  checking the results  ###################

# b need to be negative
hist(results$b)
mean(results$b)

# need to see a high positive correlation here
plot(beta[,1],colMeans(results$beta))
cor(beta[,1],colMeans(results$beta))

# results$pip should be high for primary instances (pip=1)
library(vioplot)
vioplot(results$pip[pip==F],results$pip[pip==T])

# FDRs
results$FDRs
plot(thetas,results$FDRs)
