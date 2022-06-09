#' plot_sankey
#'
#' @param Pathway_betas_path a character value representing the path to the file of pathway coefficients
#' @param receiver_topN an integeter value representing the number of top receiving pathways to plot
#' @param sender_topN an integeter value representing the number of top sending pathways to plot
#' @param pathway_shorten a boolean value representing whether to shorten the pathway name
#' 
#' @export
#'
#' @examples
#' 
#' Pathway_betas_path = "/Users/shijiazhu/Documents/MyPackages/Git/spacia/data_sz/shijia_simulation/Pathway_betas.csv"
#' plot_sankey(Pathway_betas_path)
#' plot_sankey(Pathway_betas_path, receiver_topN=5, sender_topN=3, pathway_shorten=FALSE)
#' 
plot_sankey <- function(Pathway_betas_path="Pathway_betas.csv", 
                        receiver_topN=5, sender_topN=3, 
                        pathway_shorten=TRUE) {
    
    ######## reorganize the input file
    Pathway_betas = read.csv(Pathway_betas_path, header=T)
    receiver = sapply(strsplit(as.character(Pathway_betas[,1]), "_"), 
                      function(x) paste(x[1],x[2],sep="."))
    sender = as.character(Pathway_betas[,2])    
    beta = Pathway_betas[,3]  
    data = data.frame(receiver, sender, beta)
    
    ######## scale and choose topN
    sumR = with(data, tapply(beta, receiver, function(x) sum(abs(x)) ))
    indR = names(sumR)[order(sumR, decreasing=T)][1:min(receiver_topN,length(sumR))]
    dataR = subset(data, receiver%in%indR)
    dataSR = do.call(rbind, lapply(split(dataR, as.character(dataR$receiver)), function(x) {
        y = data.frame( x, beta_scaled = abs(x$beta)/sum(abs(x$beta)) )
        y[order(y$beta_scaled, decreasing=T), ][1:min(sender_topN,nrow(y)), ]
    }))
    
    ######## get unique pathways
    if(pathway_shorten) {
        
        uni_rp = unique(as.character(dataSR$receiver))
        uni_sp = unique(as.character(dataSR$sender))
        rps = data.frame( original=uni_rp, replace=paste0("rp",1:length(uni_rp)), stringsAsFactors=F )
        sps = data.frame( original=uni_sp, replace=paste0("sp",1:length(uni_sp)), stringsAsFactors=F )
        
        p = function(ps) print(paste( paste(ps[,2],ps[,1],sep="->"),collapse="; "))
        cat("Receiving pathways:\n")
        p(rps)
        cat("Sending pathways:\n")
        p(sps)
        
        dataSR = dataSR
        dataSR$sender = with( sps, replace[ match( as.character(dataSR$sender), original ) ] )
        dataSR$receiver = with( rps, replace[ match( as.character(dataSR$receiver), original ) ] )
    }
    
    plot_interaction_river(dataSR, Bname="sender", Aname="receiver", Wname ="beta_scaled")
    
}



plot_interaction_river <- function(data, 
                       Aname=colnames(data)[1], Bname=colnames(data)[2], Wname=colnames(data)[3], 
                       Acolor = NULL, Bcolor = NULL, 
                       title = NULL ) {
    library(ggplot2)
    library(ggalluvial)
    
    data = data.frame( as.character(data[,Aname]),
                       as.character(data[,Bname]),
                       as.numeric(data[,Wname]), 
                       stringsAsFactors = F)
    colnames(data) = c(Aname,Bname,Wname)
    
    scPalette <- function(n) 
    {
        colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                        "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
                        "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
                        "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                        "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
                        "#8DD3C7", "#999999")
        if (n <= length(colorSpace)) {
            colors <- colorSpace[1:n]
        }
        else {
            colors <- (grDevices::colorRampPalette(colorSpace))(n)
        }
        return(colors)
    }
    
    ggPalette <- function(n) 
    {
        hues = seq(15, 375, length = n + 1)
        grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    ABname <- c(Aname, Bname)
    
    axes = match(ABname, colnames(data))
    data.lodes <- to_lodes_form(data, axes = axes, id = "connection")

    nA =  length(unique( data[, Aname] ))
    nB =  length(unique( data[, Bname] ))
    
    if (is.null(Acolor)) {
        Acolor <- scPalette(nA)
    }
    if (is.null(Bcolor)) {
        Bcolor <- ggPalette(nB)
    }
    ABcolor <- c(Acolor, Bcolor)
    
    font.size = 2.5
    font.size.title = 12
    col.direction = "forward"
    
    gg <- ggplot( data.lodes, aes(x = factor(x, levels = ABname), y = get(Wname), 
                                  stratum = stratum, alluvium = connection, fill = stratum, label = stratum)) + 
        geom_flow(width = 0.02, aes.flow = col.direction) + 
        geom_stratum(width = 0.02, size = 0.1, color = "black", alpha = 0.8, linetype = 1) + 
        geom_text(stat = "stratum", aes(label = after_stat(stratum)),hjust=-0.3) + #hjust: move 'brick' labels to right
        scale_x_discrete(limits = c(), labels = c("Receiving pathways", "Sending pathways")) + 
        scale_fill_manual(values = alpha(ABcolor, alpha = 0.8), drop = FALSE) + 
        theme_bw() + 
        theme(legend.position = "none", axis.title = element_blank(), axis.text.y = element_blank(), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 10)) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        ggtitle(title)

    
    gg
    
}




