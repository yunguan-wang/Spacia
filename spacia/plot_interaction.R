#' plot_interaction
#'
#' @param sending_pathway the pathway of interest to plot. If sending_pathway=="all", then draw all pathways by aggregating the signals.
#' @param receiving_pathway_path a character value representing the path to receiving pathway result folder
#'
#' @return
#' @export
#'
#' @examples
#' receiving_pathway_path = "/Users/shijiazhu/Documents/MyPackages/Git/spacia/data/simulation"
#' plot_interaction(receiving_pathway_path, sending_pathway="all")
#' plot_interaction(receiving_pathway_path, sending_pathway="1")
#' 
plot_interaction <- function(receiving_pathway_path,
                             sending_pathway) {   
    
    library(rjson)

    my_col2rgb = function(col, transparency)
    {
        color = col2rgb(col)
        rgb(color[1],color[2],color[3],max=255,alpha=transparency*255)
    }
    
    prior_path = setwd(receiving_pathway_path)
    metadata_path = "simulation_metadata.txt"
    exp_receiver_path = "exp_receiver.csv"
    exp_sender_path = "exp_sender.json"
    beta_path = "_beta.txt"
    sender_col='brown2'
    receiver_col='dodgerblue1'
    interaction_col='springgreen2'
    PI_size = 1
    interaction_size=0
    
    ########### read spot location
    loc = read.delim(metadata_path, sep="\t", header=T, row.names=1)
    
    all_sender = rownames(subset(loc,Celltype=="Senders"))
    all_receiver = rownames(subset(loc,Celltype=="Receivers"))
    
    PI_inter = subset(loc, Sender_cells_PI!="")
    PI_sender = strsplit(as.character(PI_inter$Sender_cells_PI),",")
    PI_receiver = rownames(PI_inter)
    
    ########### read gene expression
    receivers_mtx = read.csv(exp_receiver_path, header=F, row.names = NULL, stringsAsFactors = F)
    exp_receiver = receivers_mtx$V1 
    exp_sender = fromJSON(file=exp_sender_path)
    exp_sender = sapply(exp_sender, function (x) do.call(rbind, x))
    
    ########### find the expr matrix
    tmp_sender = subset(loc, Sender_cells!="")
    tmp_sender = strsplit(as.character(tmp_sender$Sender_cells),",")
    exp_sender2 = do.call(rbind, exp_sender)
    index = tapply( 1:nrow(exp_sender2), do.call(c, tmp_sender), function(x) x[1] )
    exp_sender2 = exp_sender2[index,]
    rownames(exp_sender2) = names(index)
    colnames(exp_sender2) = 1:ncol(exp_sender2)
    exp_sender2 = exp_sender2[order(as.integer(rownames(exp_sender2))), ]
    
    tmp_receiver = subset(loc, Celltype=="Receivers" & Sender_cells!="")
    index = which(tmp_receiver$Sender_cells_PI!="")
    exp_receiver2 = exp_receiver[ index ]
    
    ########### read and merge beta
    beta = read.delim(beta_path, sep="\t", header=T)
    betaS = apply(beta, 2, function(x) mean(x[ x>quantile(x,0.1) & x<quantile(x,0.9) ]) ) 
    
    PIs = unique(c( PI_receiver, do.call(c,PI_sender) ))
    all_sender = setdiff(all_sender, PIs)
    all_receiver = setdiff(all_receiver, PIs)
    
    ########### plot all nodes
    Pmargin = par()$mar
    par(mar=c(0,0,0,0))
    plot.new()
    plot.window(xlim=range(loc$X), ylim=range(loc$Y))
    points(loc[, c("X","Y")], cex=0.5, pch=16, col="darkgrey" )
    points(loc[all_sender, c("X","Y")], cex=0.5, pch=16, col=sender_col )
    points(loc[all_receiver, c("X","Y")], cex=0.5, pch=16, col=receiver_col )
    
    ########### plot PI 
    for(j in 1:length(PI_receiver))
    {
        
        senderj = PI_sender[[j]]
        xj = loc[ PI_receiver[j], "X" ]
        yj = loc[ PI_receiver[j], "Y" ]
        
        sdj = c()
        for( i in 1:length(senderj))
        {
            if(sending_pathway=='all') {
                sdj[i] = sum( betaS * exp_sender2[ senderj[i], ] )
            } else {
                exp_p = exp_sender2[ , sending_pathway ]
                sdj[i] = exp_p[ senderj[i] ]
            }
        }   
            
        for( i in 1:length(senderj))
        {
            wij = pnorm( sdj[i] )
            xi = loc[ senderj[i], "X" ]
            yi = loc[ senderj[i], "Y" ]
            ci = my_col2rgb(sender_col,wij)
            arrows(xi , yi , xj, yj, col=interaction_col, length=interaction_size, lwd=1 )
            points( xi, yi, col=ci, pch=16, cex=PI_size)
        }
        
        if(sending_pathway=='all') {
            rdj = 1 - abs( exp_receiver2[j] - pnorm(sum(sdj)) )
        } else {
            rdj = exp_receiver2[j]
        }
        
        cj = my_col2rgb(receiver_col, rdj)
        points(xj,yj, col=cj, pch=16, cex=PI_size)
    }

    
    par(mar=Pmargin)
    setwd(prior_path)
}

