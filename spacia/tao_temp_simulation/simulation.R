set.seed(1)

n_cell=5000 # simulate cells
coord=as.data.frame(matrix(0,ncol=2,nrow=n_cell))
colnames(coord)=c("X","Y")
rownames(coord)=paste("cell",1:dim(coord)[1],sep="")
coord$X=runif(n_cell,0,1)
coord$Y=runif(n_cell,0,1)

coord$type="sender" # some are senders and some are receivers
coord$type[sqrt((coord$X-0.2)^2+(coord$Y-0.3)^2)<0.1]="receiver"
coord$type[sqrt((coord$X-0.8)^2+(coord$Y-0.9)^2)<0.15]="receiver"
coord$type[sqrt((coord$X-0.2)^2+(coord$Y-0.7)^2)<0.05]="receiver"
coord$type[sqrt((coord$X-0.8)^2+(coord$Y-0.2)^2)<0.2]="receiver"
coord$type[sqrt((coord$X-0.5)^2+(coord$Y-0.5)^2)<0.05]="receiver"
coord$type[sqrt((coord$X-0.1)^2+(coord$Y-0.8)^2)<0.05]="receiver"
plot(coord$X,coord$Y,cex=0.1,col=as.numeric(as.factor(coord$type))+1)

# all 50 columns (pathways) of the senders will impact the first 
# column of the receiver
n_features=50
exp_data=matrix(rnorm(n_features*dim(coord)[1]),
                ncol=n_features,nrow=dim(coord)[1])
coord$bag=""
coord$primary=""

# assume the first three pathways in senders have the strongest potential
# of impacting the first pathway in the receivers
beta=rnorm(n_features)
hist(beta)
beta[1:3]=c(5,10,-10)
exp_sender=pos_sender=list()

for (i in 1:dim(coord)[1])
{
    if (coord$type[i]=="sender") {next}
    distance=sqrt((coord$X-coord$X[i])^2+(coord$Y-coord$Y[i])^2)
    
    # define bags of instances
    keep=distance<0.05 & coord$type=="sender"
    if (sum(keep)==0) {next}
    coord$bag[i]=paste(which(keep),collapse=",")
    exp_sender[[i]]=exp_data[keep,,drop=F]
    pos_sender[[i]]=distance[keep]
    
    # define primary instances
    keep=distance<0.04 & coord$type=="sender"
    coord$primary[i]=paste(which(keep),collapse=",")
    
    # simulate receiver exp
    exp_data[i,1]=rnorm(1,0,0.1)
    if (sum(keep)==0) {next}
    exp_data[i,1]=exp_data[i,1]+sum(colSums(exp_data[keep,,drop=F])*beta)
}

# some final organizations
exp_receiver=exp_data[1:length(exp_sender),1]>0
names(exp_sender)=names(pos_sender)=
    paste("cell",1:length(exp_sender),sep="")
keep=!sapply(pos_sender,function(x) is.null(x))
exp_receiver=exp_receiver[keep]
exp_sender=exp_sender[keep]
pos_sender=pos_sender[keep]

length(exp_receiver)
hist(sapply(strsplit(coord[names(exp_sender),"bag"],","),
            function(x) length(x)),main="bag size distribution")
hist(sapply(strsplit(coord[names(exp_sender),"primary"],","),
            function(x) length(x)),main="#primary instance distribution")

save(exp_receiver,exp_sender,pos_sender,
     beta,coord,file="~/iproject/simulated.RData")
