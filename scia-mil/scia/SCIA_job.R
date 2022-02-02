# This is a wrapper script for MIL_wrapper.R, serving as the interface between
# python master script and the MIL model.

library(Rcpp)

# Debug data
setwd('E:/projects/cell2cell_inter/data/')
scia_path = '../code/scia-mil/scia/'
job_cts = 'job_cts.txt'
job_dist = 'job_dist.txt'
job_receiver = 'job_receiver.txt'

######## Setting up ########

args = commandArgs(trailingOnly=TRUE)
scia_path = args[1]
job_cts = args[2]
job_dist = args[3]
job_receiver = args[4]

#########  source codes  #################
sourceCpp(paste(scia_path,"Fun_MICProB_C2Cinter.cpp", sep=''))
source(paste(scia_path,'MICProB_MIL_C2Cinter.R', sep=''))
source(paste(scia_path,'MIL_wrapper.R', sep=''))

######## format input into proper formats ########
job_cts = read.table(job_cts)
job_dist = read.table(job_dist, header=T, row.names = 1)
job_receiver = read.table(job_receiver, header=T, row.names = 1)

exp_receiver = job_receiver$exp_receiver

