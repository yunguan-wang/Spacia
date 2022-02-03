# This is a wrapper script for MIL_wrapper.R, serving as the interface between
# python master script and the MIL model. The only usage for this script so far 
# is to be called by the python wrapper.

library(Rcpp)

# Debug data
# setwd('E:/projects/cell2cell_inter/data/')
# scia_path = '../code/scia-mil/scia/'
# job_cts = 'job_cts.txt'
# job_dist = 'job_dist.txt'
# job_receiver = 'job_receiver.txt'

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
scia_path = args[1]
job_cts = args[2]
job_dist = args[3]
job_receiver = args[4]
job_id = args[5]
ntotal = as.integer(args[6])
nwarm= as.integer(args[7])
nthin= as.integer(args[8])
nchain= as.integer(args[9])
output_path = args[10] # output path need to have '/' at the end
thetas = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
# redirect logs
sink(file = paste(output_path, job_id, '_log.txt', sep=''))
#########  source codes  #################
sourceCpp(paste(scia_path,"Fun_MICProB_C2Cinter.cpp", sep=''))
source(paste(scia_path,'MICProB_MIL_C2Cinter.R', sep=''))
source(paste(scia_path,'MIL_wrapper.R', sep=''))

######## format input into proper formats ########
job_cts = read.table(job_cts,stringsAsFactors = F)
job_dist = read.table(job_dist, header=T, row.names = 1,check.names = F,stringsAsFactors = F)
job_receiver = read.table(job_receiver, header=T, row.names = 1,stringsAsFactors = F)

# Construct the receiver vector
exp_receiver = job_receiver$exp_receiver

# construct the sender positions
job_receiver_list = as.list(job_receiver)$senders
pos_sender = sapply(
  1:dim(job_receiver)[1],
  function (i) unname(
    unlist(
      job_dist[i,strsplit(job_receiver[i,1],',')[[1]]])
  ))

# construct the sender expression
exp_sender = lapply(
  job_receiver_list,
  function (x) unname(data.matrix((job_cts[strsplit(x,',')[[1]],]), rownames.force = F)))

res = MIL_C2Cinter(exp_receiver, pos_sender, exp_sender, ntotal, nwarm, nthin, nchain, thetas)

# save job result to disk
for (n in names(res)) {
    if (n == 'FDRs') {
        fdr = res$FDRs
        fdr[is.na(fdr)] = 1
        write.table(fdr, paste(output_path, job_id,'_',n,'.txt', sep=''), sep='\t')
    } else {
        write.table(res[n], paste(output_path, job_id,'_',n,'.txt', sep=''), sep='\t')
    }
}
