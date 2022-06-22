# This is a wrapper script for MIL_wrapper.R, serving as the interface between
# python master script and the MIL model. The only usage for this script so far 
# is to be called by the python wrapper.

library(Rcpp)
library(rjson)
# Debug data
# setwd('E:/projects/cell2cell_inter/code/data/simulation/')
# spacia_path = './spacia/'
# exp_receiver = 'exp_receiver.csv'
# exp_sender = 'exp_sender.json'
# dist_sender = 'dist_sender.json'
# ntotal=8000
# nwarm=4000
# nthin=10
# nchain=4
# thetas=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
# output_path = ''
# job_id = ''

# Dummy test data
# exp_receiver=c(1,1,0,0)
# pos_sender=list(
#  c(1,4,3),
#  c(1),
#  c(1,5,2,3,6),
#  c(3)
# )
# 
# # instance features, expression of genes/pathways in the sending cells
# # rows are genes/pathways, cols are instances
# exp_sender=list(
#  matrix(runif(3*3),ncol=3),
#  matrix(runif(1*3),ncol=3),
#  matrix(runif(5*3),ncol=3),
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

######## Setting up ########

args = commandArgs(trailingOnly=TRUE)
spacia_path = args[1]
exp_sender = args[2]
dist_sender = args[3]
exp_receiver = args[4]
job_id = args[5]
ntotal = as.integer(args[6])
nwarm = as.integer(args[7])
nthin = as.integer(args[8])
nchain = as.integer(args[9])
output_path = args[10] # output path need to have '/' at the end
plot_mcmc = as.logical(args[11]) # whether or not to plot diagnosis plots
ext = args[12] # extension used for the ggsave. see BetaB2MCMCPlots.R
thetas = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)


# redirect logs
sink(file = paste(output_path, job_id, '_log.txt', sep=''))
#########  source codes  #################
sourceCpp(paste(spacia_path,"Fun_MICProB_C2Cinter.cpp", sep=''))
source(paste(spacia_path,'MICProB_MIL_C2Cinter.R', sep=''))
source(paste(spacia_path,'MIL_wrapper.R', sep=''))
source(paste(spacia_path,'BetaB2MCMCPlots.R', sep=''))

######## format input into proper formats ########

# Read receiver matrix
exp_receiver = read.csv(
    exp_receiver, header=F, row.names = NULL, stringsAsFactors = F)$V1
exp_receiver = exp_receiver == 1

# Read sender expression 
exp_sender = fromJSON(file=exp_sender)
exp_sender = sapply(exp_sender, function (x) do.call(rbind, x))

# Read sender distance to receivers
dist_sender = fromJSON(file=dist_sender)
# Normalize distance with the maximal distance
max_dist = max(sapply(dist_sender, function(x) x[which.max(abs(x))]))
dist_sender = sapply(dist_sender, function(x) x / max_dist)

# Run the model 
set.seed(0)
t0 = Sys.time()
res = MIL_C2Cinter(
  exp_receiver, dist_sender, exp_sender, 
  ntotal, nwarm, nthin, nchain, thetas)
t1 = Sys.time()
print(t1-t0)
# Get memory use
gc()
# save job result to disk
for (n in names(res)) {
    if (n == 'FDRs') {
        fdr = res$FDRs
        fdr[is.na(fdr)] = 1
        write.table(
          fdr, paste(output_path, job_id,'_',n,'.txt', sep=''), sep='\t')
    } else {
        write.table(
          res[n], paste(output_path, job_id,'_',n,'.txt', sep=''), sep='\t')
    }
}

########### Plot MCMC Diagnostics ##############
if (plot_mcmc) {
  
  beta_matrix = as.matrix(res$beta)
  b_matrix = as.matrix(res$b)
  colnames(beta_matrix) = paste("beta.", 1:dim(beta_matrix)[2], sep="")
  colnames(b_matrix) = c("b.1", "b.2")
  
  S <- BetaB2MCMCPlots(beta_matrix,
                       b_matrix,
                       nwarm,
                       ntotal,
                       nthin,
                       nchain,
                       job_id,
                       output_path,
                       ext)
  
}