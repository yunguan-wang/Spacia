##########  set Rstudio environment  ############

# this is a problem only specific to the Rstudio of biohpc
Sys.setenv(LIBRARY_PATH = "/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mkl/lib/intel64:/cm/shared/apps/java/oracle/jdk1.7.0_51/lib:/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/compiler/lib/intel64:/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mpi/intel64/lib:/cm/shared/apps/openmpi/gcc/64/2.1.5/lib64:/cm/shared/apps/gcc/5.4.0/lib:/cm/shared/apps/gcc/5.4.0/lib64:/cm/shared/apps/slurm/16.05.8/lib64/slurm:/cm/shared/apps/slurm/16.05.8/lib64")

########  load simulated data  ###########

load("~/iproject/simulated.RData")

################  running the MIL model  ###################

# other input parameters
ntotal=8000
nwarm=4000
nthin=10
nchain=4

thetas=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# run MIL
# setwd("/project/shared/xiao_wang/projects/cell2cell_inter/code/spacia/spacia")
setwd("E:/projects/cell2cell_inter/code/spacia/spacia")
library('Rcpp')
source('MIL_wrapper.R')
sourceCpp("Fun_MICProB_C2Cinter.cpp")
source('MICProB_MIL_C2Cinter.R')

Sys.time()
set.seed(0)
results=MIL_C2Cinter(exp_receiver,pos_sender,exp_sender,
                     ntotal,nwarm,nthin,nchain,thetas)
Sys.time()

###########  checking the results  ###################

# b need to be negative
hist(results$b[,2])
mean(results$b[,2])

# need to see a high positive correlation here
plot(beta,colMeans(results$beta))
cor(beta,colMeans(results$beta))

pi=c()
coord_bag=coord[names(exp_sender),]

for (i in 1:dim(coord_bag)[1])
{
    tmp=strsplit(coord_bag$bag[i],",")[[1]] %in% 
        strsplit(coord_bag$primary[i],",")[[1]]
    pi=c(pi,tmp)
}

boxplot(results$pip_recal~pi,names=c("non-primary","primary"),
        xlab="",ylab="primary-instance score")

# results$pip should be high for primary instances (pip=1)
# the coding is a bit complicated here. I will figure it out tomorrow
#library(vioplot)
#vioplot(results$pip[pip==F],results$pip[pip==T])

# FDRs
results$FDRs
plot(thetas,results$FDRs)

plot(1:1604, results$beta[,5], type = 'l')
