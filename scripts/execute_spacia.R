library(Rcpp)
library(rjson)
library("optparse")

option_list = list(
  make_option(c("-i", "--inputDir"), type="character", 
              default=NULL, 
              help='input data directory [default = %default]\nmust contain:\n\tspacia_spot_meta.txt: table of location and cell type info.\n\t\tMust contain cell names and row names as well as columns "X","Y", and "cell_type"\n\tcell_by_gene.csv: gene expression matrix of cells. Columns shouldbe "cell" followed by gene names',
              metavar="character"),
  make_option(c("-a", "--spacia_path"), type="character", default=NULL, 
              help="path to spacia core code [default = %default]\nThe path to directory containing:\n\tFun_MICProB_C2Cinter.cpp\n\tMICProB_MIL_C2Cinter.R\n\tMIL_wrapper.R", 
              metavar="character"),
  make_option(c("-r", "--receivingCell"), type="character", default=NULL, 
              help="receiving cell type [default = %default]", metavar="character"),
  make_option(c("-s", "--sendingCell"), type="character", default=NULL, 
              help="sending cell type [default = %default]", metavar="character"),
  make_option(c("-g", "--receivingGene"), type="character", default=NULL, 
              help="receiving gene [default = %default]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name prefix", metavar="character"),
  make_option(c("-f", "--overwrite"), action="store_true", default=FALSE,
              help="overwrite existing output [default = %default]"),
  make_option(c("-t", "--paramTable"), type="character", 
              help="csv of cor_cutoff and exp_receiver_quantile for each receiving gene\nmust contain columns:\n\tgene_name,cor_cutoff,quantile_cutoff", 
              metavar="character"),
  make_option(c('-q', '--quantile'), type='double', default=NULL,
              help='receiving gene quantile cutoff, overwrites -t [default = %default]\nquantile_cutoff is the cutoff used to dichotomize the receiving gene signature', 
              metavar = 'number'),
  make_option(c('-u', '--corCut'), type='double', default=NULL,
              help='receiving gene cor. cutoff, overwrites -t [default = %default]\ncor_cutoff is the correlation value cutoff used in choosing genesto construct a signature of the receiving gene; this reduces dropout', 
              metavar = 'number'),
  make_option(c('-d', '--dist'), type='integer', default=50,
              help='distance cutoff [default = %default]', metavar = 'number'),
  make_option(c('-p', '--path'), type='integer', default=50,
              help='number of principle components to use [default = %default]', metavar = 'number'),
  make_option(c('-m', '--min'), type='integer', default=3,
              help='min number of instances per bag [default = %default]', metavar = 'number'),
  make_option(c('-b', '--subSample'), type='integer', default=5000,
              help='maximum number of bags [default = %default]', metavar = 'number'),
  make_option(c('-l', '--ntotal'), type='integer', default=50000,
              help='ntotal [default = %default]', metavar = 'number'),
  make_option(c('-w', '--nwarm'), type='integer', default=25000,
              help='nwarm [default = %default]', metavar = 'number'),
  make_option(c('-n', '--nthin'), type='integer', default=10,
              help='nthin [default = %default]', metavar = 'number'),
  make_option(c('-c', '--nchain'), type='integer', default=3,
              help='nchain [default = %default]', metavar = 'number'),
  make_option(c('-e', '--nSample'), type='integer', default=50,
              help='number of samples from each chain to calculate beta/b pvals [default = %default]', metavar = 'number')
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
########  parameters  ###########
cat('################starting run...################\n')
Sys.time()
if (is.null(opt$output)) {
  outFn = paste(opt$sendingCell, '-', opt$receivingCell, '_', opt$receivingGene)
} else{
  outFn = opt$output
}
rDataFn = paste(outFn, '.RData', sep = '')
if (file.exists(rDataFn)) {
  if (opt$overwrite) {
    cat('ignoring existing output...\n')
  } else{
    Sys.time()
    stop('*********terminating run: found existing output*********\n')
  }
}
# choose dataset, cell type pairs and genes to investigate
setwd(opt$inputDir)
sending_cell_type = opt$sendingCell
cat(paste('sending cells:', sending_cell_type, '\n'))
receiving_cell_type = opt$receivingCell
cat(paste('receiving cells:', receiving_cell_type, '\n'))
receiving_gene=opt$receivingGene
cat(paste('receiving gene:', receiving_gene, '\n'))
nSample = opt$nSample


# important tuning parameters, users may need to play with them
# runtime is very good with current selection of parameters
if (is.null(opt$quantile) & is.null(opt$corCut)) {
  
  paramTable = read.csv(opt$paramTable)
  tmp = paramTable[paramTable$gene_name == receiving_gene, ]
  if (dim(tmp)[1] == 1) {
    cor_cutoff = tmp$cor_cutoff
    exp_receiver_quantile = tmp$quantile_cutoff
    writeLines(paste('found ', receiving_gene, ':\n\tcor_cutoff: ', cor_cutoff, 
                     '\n\texp_receiver_quantile: ', exp_receiver_quantile, sep = ''))
  } else if (dim(tmp)[1] == 0){
    stop(paste('error:', receiving_gene, 'not found in table'))
  }else if (dim(tmp)[1] > 1){
    stop(paste('error: multiple matches for', receiving_gene, 'found in table; check table'))
  }
} else{
  cor_cutoff = opt$corCut
  cat(paste('cor cutoff:', cor_cutoff, '\n'))
  exp_receiver_quantile = opt$quantile
  cat(paste('quantile cutoff:', exp_receiver_quantile, '\n'))
}

dist_cutoff=opt$dist
cat(paste('distance cutoff:', dist_cutoff, '\n'))
n_path=opt$path
cat(paste('num. PCs:', n_path, '\n'))
min_instance=opt$min
cat(paste('min. instance/bag:', min_instance, '\n'))

# other less important input parameters
ntotal=opt$ntotal
cat(paste('ntotal:', ntotal, '\n'))
nwarm=opt$nwarm
cat(paste('nwarm:', nwarm, '\n'))
nthin=opt$nthin
cat(paste('nthin:', nthin, '\n'))
nchain=opt$nchain
cat(paste('nchain:', nchain, '\n'))
thetas=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

########  read and pre-process data  ############
# many of these steps are assumed to be done by users
# we mostly use user-provided input as is
Sys.time()
cat('reading data...\n')
options(scipen=999) 
# do this in spacia
meta=read.table(file.path(opt$inputDir, "spacia_spot_meta.txt"),stringsAsFactors = F,sep="\t")

counts=read.csv(file.path(opt$inputDir, "cell_by_gene.csv"),stringsAsFactors=F)
Sys.time()
cat('pre-processing data...\n')
#rownames(counts)=paste("cell_",counts$cell+1,sep="") # user's responsibility to name cells consistently and in a good manner
counts=as.matrix(counts[,-1])
counts=t(t(counts)/colMeans(counts)) # user's responsibility, don't do it in spacia
counts=log(counts+1) # user's responsibility, don't do it in spacia
counts=t(counts) 
counts=counts[apply(counts,1,sd)>0,] # make sure the scale of each gene's expression is the same, do it in spacia
counts=counts/apply(counts,1,sd) # make sure the scale of each gene's expression is the same, do it in spacia
counts=counts[!grepl("Blank",rownames(counts)),] # user's responsibility, don't do it in spacia

if (sum(!rownames(meta) %in% colnames(counts))>0) 
  {stop("Cell name mismatch!")} 
counts=counts[,rownames(meta)]

meta_tumor=meta[meta$cell_type==receiving_cell_type,]
counts_tumor=counts[,rownames(meta_tumor)]
cat(paste('\treceiving cells:', dim(meta_tumor)[1], '\n'))
meta_fibro=meta[meta$cell_type==sending_cell_type,]
counts_fibro=counts[,rownames(meta_fibro)]
cat(paste('\tsending cells:', dim(meta_fibro)[1], '\n'))

# we have hundreds to thousands of genes for senders 
# (fibro here, I am just too lazy to change the naming of variables)
# aggregate them to a few dozen principal components
# put through spacia
# when we get the spacia results, we will revert back to gene level
Sys.time()
cat('calculating pca...\n')
pca_fibro=prcomp(t(counts_fibro))$x
Sys.time()
cat('processing pca...\n')
pca_fibro=t(t(pca_fibro)/apply(pca_fibro,2,sd)) # do this in spacia
Sys.time()
cat('constructing bags...\n')

pos_sender=exp_sender=list()
xy_tumor=t(as.matrix(meta_tumor[,c("X","Y")]))
xy_fibro=t(as.matrix(meta_fibro[,c("X","Y")]))
Sys.time()

#find pos and exp of senders within dist_cutoff
for (i in 1:dim(counts_tumor)[2])
{
  #slow with large num. of sender cells
  dist_tumor_fibro=colSums((xy_tumor[,i]-xy_fibro)^2)^0.5
  keep=dist_tumor_fibro<dist_cutoff
  if (sum(keep)<min_instance) {next}
  pos_sender[[rownames(meta_tumor)[i]]]=log(dist_tumor_fibro[keep]) # critical (log)
  exp_sender[[rownames(meta_tumor)[i]]]=pca_fibro[keep,1:n_path]
}
nbags = length(pos_sender)
cat("Total number of receiving cells:",dim(meta_tumor)[1],"\n")
cat("Successfully constructed bags:",nbags,"\n")
Sys.time()
cat('finalizing spacia inputs...\n')
# hist(sapply(pos_sender,function(x) length(x)),main="number of instances per bag")

# we will iterate through each gene in the receivers or the genes that 
# the users specify (maybe a few dozens)
# right now spacia runs rather fast. we can handle 500 receiving
# genes if we really want to
cors=cor(counts_tumor[receiving_gene,],t(counts_tumor))[1,]
keep=abs(cors)>cor_cutoff
cat(paste(sum(keep),"genes highly correlated with",receiving_gene,"\n"))
# aggregate receiving genes to receiving pathways
signature=colMeans(counts_tumor[keep,names(exp_sender)]*cors[keep])
# hist(signature,breaks=30,main="Distribution of the expression of the signature")
# abline(v=quantile(signature,exp_receiver_quantile))
exp_receiver=1*(signature>quantile(signature,exp_receiver_quantile))
# users may want to do some trial and errors to choose the tuning parameters 
# for a good "exp_receiver" 
# for our pan-cancer spacia analyses, we may not have the luxury 
# to choose parameters for all spacia runs 
# across tumor cells, gene pairs, etc. 
# for individual analyses, definitely tune
maxBags = opt$subSample
if (maxBags > 0) {
  if (nbags > maxBags) {
    cat("Subsampling constructed bags to",maxBags,"\n")
    keep1 = sample(1:length(exp_receiver), maxBags)
    exp_receiver = exp_receiver[keep1]
    pos_sender = pos_sender[keep1]
    exp_sender = exp_sender[keep1]
  }
}
Sys.time()



######  run spacia  #####################

cat('loading spacia...\n')
#path for required libraries for spacia (uncomment and edit if needed)
#Sys.setenv(LIBRARY_PATH = "/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mkl/lib/intel64:/cm/shared/apps/java/oracle/jdk1.7.0_51/lib:/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/compiler/lib/intel64:/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mpi/intel64/lib:/cm/shared/apps/openmpi/gcc/64/2.1.5/lib64:/cm/shared/apps/gcc/5.4.0/lib:/cm/shared/apps/gcc/5.4.0/lib64:/cm/shared/apps/slurm/16.05.8/lib64/slurm:/cm/shared/apps/slurm/16.05.8/lib64")
spacia_path = opt$spacia_path
sourceCpp(paste(spacia_path,"Fun_MICProB_C2Cinter.cpp", sep=''))
source(paste(spacia_path,'MICProB_MIL_C2Cinter.R', sep=''))
source(paste(spacia_path,'MIL_wrapper.R', sep=''))

cat('running spacia...\n')
Sys.time()
results=MIL_C2Cinter(exp_receiver,pos_sender,exp_sender,
                     ntotal,nwarm,nthin,nchain,thetas, 1)
Sys.time()

cat('processing results...\n')
# this is important
# rotation: the contribution of each sending gene expression to each PC
# beta: the contribution of each PC to receiving gene expression
# sum(rotation*beta): the contribution of each sending gene
#                     to each receiving gene
rotation=prcomp(t(counts_fibro))$rotation
# remember we have a lot of beta from multiple MCMC iterations
# don't just take mean(beta) and multiple with rotation
# multiply each beta with rotation and then sum. Then
# these sums will form a distribution from which we can do stat test
gene_level_beta=results$beta %*% t(rotation[,1:dim(results$beta)[2]])

#####  examine results  ########
# check top genes with the largest/smallest gene-level beta
tmp=colMeans(gene_level_beta)
cat('top 10 gene betas')
tmp[rank(-abs(tmp))<10]

# needs to negative
cat(paste('mean b:', mean(results$b), '\n'))

Sys.time()
cat('saving raw results...\n')

save(gene_level_beta, results, nbags, exp_receiver,pos_sender,exp_sender,
     sending_cell_type , receiving_cell_type, receiving_gene,
     file = rDataFn)
cat(paste('saved results to', rDataFn, '\n'))
cat('saving csv results...\n')

#recover order from structure of input 
#EX:
#	{receivingCell1: {sendingCell1: val, sendingCell2: val...}, 
#  receivingCell2: {sendingCell1: val, sendingCell3: val...}, ...}
#
options(scipen=999)
sendL = c()
recL = c()
for (rec in names(pos_sender)) {
  tmp = names(pos_sender[[rec]])
  recL = c(recL, rep(rec, length(tmp)))
  sendL = c(sendL, tmp)
}

#sample nSample of each MCMC chain to calculate beta/b pvals
##check results
nOutPerChain = (ntotal - nwarm) / nthin + 1
nOutExpected = nOutPerChain * nchain

nBeta = dim(gene_level_beta)[1]
if (nBeta != nOutExpected) {
  stop(paste('error:', basename(outFn), 'contains', nBeta, 'betas', nOutExpected, 'expected.'))
}
nB = dim(results$b)[1]
if (nB != nOutExpected) {
  stop(paste('error:', basename(outFn), 'contains', nB, "b's", nOutExpected, 'expected.'))
}
##get pvals
betas0 = c()
bs = c()
ind = 0
ii = c()
for (n in 1:nchain) {
  #start at ~20% of the MCMC cycles
  i = round(nOutPerChain / 5)
  # inds = sample((ind + 1): (ind + nOutPerChain), nSample)
  #evenly sample to the end of each chain
  offsetI = (n -1) * nOutPerChain
  inds = round(seq.int(i + offsetI, nOutPerChain + offsetI, length.out = 50))
  ii = c(ii, inds)
  betas0 = rbind(betas0, gene_level_beta[inds,])
  bs = c(bs,results$b[inds, 2])
}
###one sided t-test for b (expect b<0 for true interactions)
testRes = t.test(bs, alternative = 'less')
bPval = testRes$p.value
###two sided t-test for beta
testRes = apply(betas0, 2, function(x){res = t.test(x); res$p.value})
#simple avg. for b; col 2 is used for this version of spacia
# b = colMeans(results$b)[2]
# b = mean(bs)
#average beta over each iteration (row)
#betas = colMeans(gene_level_beta)
# betas = colMeans(betas0)
# bl = length(betas)
betas = data.frame('sending_gene' = names(testRes),
                   'receiving_gene' = receiving_gene,
                   'avg_beta' = colMeans(gene_level_beta),
                   'avg_beta_sampled' = colMeans(betas0),
                   "beta_pval" = testRes,
                   'b' = colMeans(results$b)[2],
                   "b_sampled" = mean(bs),
                   "b_pval" = bPval)

df = data.frame('receiving_cell' = recL,
                'sending_cell' = sendL, 
                'receiving_gene' = rep(receiving_gene, length(sendL)),
                'avg_primary_instance_score' = rowMeans(results$pip)) #each col is from one mcmc chain
save(gene_level_beta, results, nbags, exp_receiver,pos_sender,exp_sender,
     ii, betas0, bs, file = rDataFn)
write.csv(betas, file = paste(outFn, '_betas.csv', sep = ''), quote = FALSE, row.names = FALSE)
write.csv(df, file = paste(outFn, '_pip.csv', sep = ''), quote = FALSE, row.names = FALSE)
cat('saved csv results\n')
Sys.time()
cat('################finished run################\n\n\n')