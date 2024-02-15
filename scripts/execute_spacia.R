library(Rcpp)
library(rjson)
library("optparse")
library(ggplot2)
library(patchwork)
library(scales)
library(gridExtra)
library(dplyr)

option_list = list(
  make_option(c("-x", "--inputExpression"), type="character", 
              default=NULL, 
              help='gene expression matrix [default = %default]\n\t\t\tFormat: cell x gene expression values. Column names should be "cell" followed by gene names. Same as in spacia.py, expression values should be normalized and log-transformed. Option "-C" can be used if using raw counts to perform a simple transformation',
              metavar="character"),
  make_option(c("-m", "--inputMeta"), type="character", 
              default=NULL, 
              help='input metadata [default = %default]\n\t\t\tTable of location and cell type info with each cell taking one row. Must contain cell names as first column as well as columns "X","Y", and "cell_type"',
              metavar="character"),
  make_option(c("-C", "--isCount"), action="store_true", 
              help="gene expression matrix consists of raw counts"),
  make_option(c("-a", "--spacia_path"), type="character", default=NULL, 
              help="path to spacia core code [default = %default]\n\t\tThe path to directory containing:\n\t\t\tFun_MICProB_C2Cinter.cpp\n\t\t\tMICProB_MIL_C2Cinter.R\n\t\t\tMIL_wrapper.R", 
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
  make_option(c("-t", "--paramTable"), type="character", default=NULL,
              help="optional csv of cor_cutoff and exp_receiver_quantile for each receiving gene\n\t\tmust contain columns:\n\t\t\tgene_name,cor_cutoff,quantile_cutoff", 
              metavar="character"),
  make_option(c('-q', '--quantile'), type='double', default=NULL,
              help='receiving gene quantile cutoff, overwrites -t [default = %default]\n\t\tcutoff used to dichotomize the receiving gene signature', 
              metavar = 'number'),
  make_option(c('-u', '--corCut'), type='double', default=NULL,
              help='receiving gene cor. cutoff, overwrites -t [default = %default]\n\t\tcorrelation value cutoff used in choosing genes to construct a signature of the receiving gene; this reduces dropout', 
              metavar = 'number'),
  make_option(c('-d', '--dist'), type='integer', default=50,
              help='distance cutoff [default = %default]', metavar = 'number'),
  make_option(c('-p', '--path'), type='integer', default=50,
              help='number of principle components to use [default = %default]', metavar = 'number'),
  make_option(c('-i', '--min'), type='integer', default=3,
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

dir.create(dirname(outFn), showWarnings = FALSE, recursive = T)
rDataFn = paste(outFn, '.RData', sep = '')
if (file.exists(rDataFn)) {
  if (opt$overwrite) {
    cat('ignoring existing output...\n')
  } else{
    Sys.time()
    stop('*********terminating run: found existing output*********\n')
  }
}
#load cached data
loadedCache = FALSE
cacheFn = paste(outFn, '_cache.RData', sep = '')
if (file.exists(cacheFn)) {
  if (opt$overwrite) {
    cat('ignoring existing cache...\n')
  } else{
    Sys.time()
    cat('loading from existing cache...\n')
    load(cacheFn)
    loadedCache = TRUE
  }
}
# choose dataset, cell type pairs and genes to investigate
if (loadedCache) {
  tmpCheck = c(sending_cell_type != opt$sendingCell,
               receiving_cell_type != opt$receivingCell,
               receiving_gene != opt$receivingGene,
               n_path != opt$path
               )
  if (any(tmpCheck)) {
    stop('*********terminating run: input mismatch with cache file; use -f to ignore cache*********\n')
  }
} 
sending_cell_type = opt$sendingCell
cat(paste('sending cells:', sending_cell_type, '\n'))
receiving_cell_type = opt$receivingCell
cat(paste('receiving cells:', receiving_cell_type, '\n'))
receiving_gene=opt$receivingGene
cat(paste('receiving gene:', receiving_gene, '\n'))
nSample = opt$nSample
runPlotOnly = FALSE
# important tuning parameters, users may need to play with them
if (is.null(opt$quantile) & is.null(opt$corCut)) {
  if (is.null(opt$paramTable)) {
    #run plots to help determine cutoffs
    runPlotOnly = TRUE
  } else {
    #if using csv
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
if (!loadedCache) {
  ########  read and pre-process data  ############
  Sys.time()
  cat('reading data...\n')
  options(scipen=999) 
  # do this in spacia
  meta=read.csv(opt$inputMeta, row.names = 1, stringsAsFactors=F)
  counts=read.csv(opt$inputExpression, row.names = 1, stringsAsFactors=F)
  Sys.time()
  cat('pre-processing data...\n')
  #process count data
  if (opt$isCount) {
    counts=as.matrix(counts[,-1])
    counts=t(t(counts)/colMeans(counts)) # norm
    counts=log1p(counts) # log
  }
  counts=t(counts) 
  # make sure the scale of each gene's expression is the same
  counts=counts[apply(counts,1,sd)>0,] 
  counts=counts/apply(counts,1,sd) 
  #check cell names
  if (sum(!rownames(meta) %in% colnames(counts))>0) 
    {stop("Cell name mismatch!")} 
  counts=counts[,rownames(meta)]
  #filter to correct cell types
  meta_receiver=meta[meta$cell_type==receiving_cell_type,]
  counts_receiver=counts[,rownames(meta_receiver)]
  cat(paste('\treceiving cells:', dim(meta_receiver)[1], '\n'))
  meta_sender=meta[meta$cell_type==sending_cell_type,]
  counts_sender=counts[,rownames(meta_sender)]
  cat(paste('\tsending cells:', dim(meta_sender)[1], '\n'))
  
  # 1. first aggregate all potential sedning genes to a few dozen principal components
  # 2. run spacia
  # 3. revert spacia results back to gene level
  Sys.time()
  cat('calculating pca...\n')
  pca_sender=prcomp(t(counts_sender))$x
  Sys.time()
  cat('processing pca...\n')
  pca_sender=t(t(pca_sender)/apply(pca_sender,2,sd)) # do this in spacia
  Sys.time()
  cat('constructing bags...\n')
  pos_sender=exp_sender=list()
  xy_receiver=t(as.matrix(meta_receiver[,c("X","Y")]))
  xy_sender=t(as.matrix(meta_sender[,c("X","Y")]))
  Sys.time()
  #find pos and exp of senders within dist_cutoff
  dist_cutoff2 = dist_cutoff^2
  for (i in 1:dim(counts_receiver)[2])
  {
    dist_receiver_sender=colSums((xy_receiver[,i]-xy_sender)^2)
    keep=dist_receiver_sender<dist_cutoff2
    if (sum(keep)<min_instance) {next}
    pos_sender[[rownames(meta_receiver)[i]]]=log(dist_receiver_sender[keep]) # critical (log)
    exp_sender[[rownames(meta_receiver)[i]]]=pca_sender[keep,1:n_path]
  }
  nbags = length(pos_sender)
  cat("Total number of receiving cells:",dim(meta_receiver)[1],"\n")
  cat("Successfully constructed bags:",nbags,"\n")
  Sys.time()
}

if (runPlotOnly) {
  save(counts_receiver,counts_sender, pos_sender, exp_sender, pca_sender, 
       nbags, sending_cell_type, receiving_cell_type, receiving_gene, 
       ntotal, nwarm, nthin, nchain, thetas, n_path,
       file = cacheFn)
  plotCutoffs <- function(receiving_gene, file_path, counts_receiver, 
                          exp_sender, pca_sender, n_path) {
    quantile_cutoffs = 1:12
    quantile_cutoffs = (quantile_cutoffs - 1) / 11
    quantile_cutoffs = quantile_cutoffs[2:11]
    # pca_cum = summary(pca_sender)$importance[3,n_path]*100
    df <- data.frame(matrix(ncol = 4, nrow = 0))
    x <- c("signature", "cor_cutoff", "receiver_quantile", "xintercept")
    colnames(df) <- x
    cors=cor(counts_receiver[receiving_gene,],t(counts_receiver))[1,]
    # Find the combination that yields 5~20 highly correlated genes
    searchgrid = c(1:1000)/1000
    min_cor_cutoff = 1
    max_cor_cutoff = 0
    for (i in searchgrid){
      num_cor_genes = sum(abs(cors) > i)
      if (num_cor_genes >= 2 & 
          num_cor_genes <= 20 ){
        if (i < min_cor_cutoff){
          min_cor_cutoff = i
        }
        if (i > max_cor_cutoff){
          max_cor_cutoff = i
        }
      }
    }
    cor_cutoffs = seq(from = min_cor_cutoff, to = max_cor_cutoff, length.out = 10)
    ncgvec = c()
    for (i in 1:10){ncgvec = c(ncgvec,sum(abs(cors)>cor_cutoffs[i]))}
    cor_df = data.frame(cor_cutoffs = cor_cutoffs,
                        num_cor_genes = ncgvec)
    for (i in cor_cutoffs){
      keep=abs(cors)>i
      signature=colMeans(counts_receiver[keep,names(exp_sender)]*cors[keep])
      #exp_receiver=sum(1*(signature>quantile(signature,j)))
      for (j in quantile_cutoffs){
        df_temp = data.frame(signature = signature, 
                             cor_cutoff = rep(i, length(signature)),
                             receiver_quantile = rep(j, length(signature)),
                             xintercept = rep(quantile(signature,j), length(signature))
        )
        df = bind_rows(df, df_temp)
      }
    }
    is.num <- sapply(df, is.numeric)
    df[is.num] <- lapply(df[is.num], round, 5)
    # Histogram with density plot
    p = ggplot(df, aes(x=signature)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 30)+
      geom_density()+
      theme_bw() +
      geom_vline(data = df, aes(xintercept = xintercept), colour = "red")+
      ggtitle(
        paste(sending_cell_type,
              " to ",
              receiving_cell_type,
              "; ",
              receiving_gene,
              "; ",n_path,"PCs",
              # pca_cum,
              sep ="")
      ) +
      facet_grid(cor_cutoff~receiver_quantile, labeller = label_both) 
    combined_plot <-p + gridExtra::tableGrob(cor_df)+ plot_layout(widths = c(5, 1))
    ggsave(file_path, plot = combined_plot, width = 40, height = 20, dpi = 500)
  }
  file_path = paste(outFn, '_cutoffPlot', '.pdf', sep = '')
  plotCutoffs(receiving_gene, file_path, counts_receiver, 
              exp_sender, pca_sender, n_path)
  s = paste('finished plotting cutoff plots to ', file_path, '')
  stop(s)
}

cat('finalizing spacia inputs...\n')
# aggregate receiving genes to receiving pathways
cors=cor(counts_receiver[receiving_gene,],t(counts_receiver))[1,]
keep=abs(cors)>cor_cutoff
cat(paste(sum(keep, na.rm = T),"genes highly correlated with",receiving_gene,"\n"))
signature=colMeans(counts_receiver[keep,names(exp_sender)]*cors[keep])
exp_receiver=1*(signature>quantile(signature,exp_receiver_quantile,na.rm = T))
# users may want to do some trial and errors to choose the tuning parameters 
# for a good "exp_receiver" 
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
#path for required libraries for spacia (edit if needed)
#Sys.setenv(LIBRARY_PATH = "/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mkl/lib/intel64:/cm/shared/apps/java/oracle/jdk1.7.0_51/lib:/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/compiler/lib/intel64:/cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mpi/intel64/lib:/cm/shared/apps/openmpi/gcc/64/2.1.5/lib64:/cm/shared/apps/gcc/5.4.0/lib:/cm/shared/apps/gcc/5.4.0/lib64:/cm/shared/apps/slurm/16.05.8/lib64/slurm:/cm/shared/apps/slurm/16.05.8/lib64")
spacia_path = opt$spacia_path
sourceCpp(file.path(spacia_path,"Fun_MICProB_C2Cinter.cpp"))
source(file.path(spacia_path,'MICProB_MIL_C2Cinter.R'))
source(file.path(spacia_path,'MIL_wrapper.R'))

cat('running spacia...\n')
Sys.time()
results=MIL_C2Cinter(exp_receiver,pos_sender,exp_sender,
                     ntotal,nwarm,nthin,nchain,thetas, 1)
Sys.time()

cat('processing results...\n')
# important:
# rotation: the contribution of each sending gene expression to each PC
# beta: the contribution of each PC to receiving gene expression
# sum(rotation*beta): the contribution of each sending gene to each receiving gene
rotation=prcomp(t(counts_sender))$rotation
# we have many beta values from multiple MCMC iterations
# instead of using mean(beta) and multiple with rotation
# we multiply each beta with rotation and then sum. Then
# these sums will form a distribution from which we can do stat test
gene_level_beta=results$beta %*% t(rotation[,1:dim(results$beta)[2]])

#####  examine results  ########
# check top genes with the largest/smallest gene-level beta
tmp=colMeans(gene_level_beta)
cat('top 10 gene betas')
tmp[rank(-abs(tmp))<10]

# needs to be negative to make sense
cat(paste('mean b:', mean(results$b), '\n'))

Sys.time()
cat('saving raw results...\n')

save(gene_level_beta, results, nbags, exp_receiver,pos_sender,exp_sender,
     sending_cell_type , receiving_cell_type, receiving_gene,
     file = rDataFn)
cat(paste('saved results to', rDataFn, '\n'))
cat('saving csv results...\n')

######  process spacia outputs  #####################
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
  i = 1
  offsetI = (n -1) * nOutPerChain
  inds = round(seq.int(i + offsetI, nOutPerChain + offsetI, length.out = 50))
  ii = c(ii, inds)
  betas0 = rbind(betas0, gene_level_beta[inds,])
  bs = c(bs,results$b[inds, 2])
}
###one sided t-test for b (expect b < 0 for true interactions)
testRes = t.test(bs, alternative = 'less')
bPval = testRes$p.value
###two sided t-test for beta
testRes = apply(betas0, 2, function(x){res = t.test(x); res$p.value})
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
if (file.exists(cacheFn)) {
  unlink(cacheFn)
}
cat('################finished run################\n\n\n')
