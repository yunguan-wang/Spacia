#' plot_b
#'
#' @param b_path a character value representing the path to the file of b
#' @param theta_cutoff a numeric value representing the cutoff of theta
#' 
#' @export
#'
#' @examples
#' 
#' setwd("/Users/shijiazhu/Documents/MyPackages/Git/spacia/data_sz/shijia_simulation/")
#' b_path = "/Users/shijiazhu/Documents/MyPackages/Git/spacia/data_sz/shijia_simulation/B_and_FDR.csv"
#' plot_b(b_path)
#' 
plot_b <- function(b_path="B_and_FDR.csv", theta_cutoff=0.9) {
    
    library(ggplot2)
    library(ggalluvial)
    library(plotrix)
    
    my_col2rgb = function(col, transparency)
    {
        color = col2rgb(col)
        rgb(color[1],color[2],color[3],max=255,alpha=transparency*255)
    }
    
    B = read.csv(b_path, header=T)
    B = B[ B$Theta_cutoff==theta_cutoff , ]
    n = nrow(B)
    
    if(1) {
        B$b = rnorm(n,sd=1.5)
        B$FDR = 10^(-abs(B$b))
    }
    
    receiver = sapply(strsplit(as.character(B[,1]), "_"), 
                      function(x) paste(x[1],x[2],sep="."))
    
    minuslog10FDR = -log10(B$FDR)
    
    # barplot
    space = 3
    xmax = (length(B$b)+2)*(space+1) - 0.5
    xlim=c(0, xmax) 
    ylim = range( B$b + sign(B$b)*minuslog10FDR)
    ylim = c( min(ylim,0), max(0,ylim) )
    col = my_col2rgb("red", 0.5 )
    
    cex = 1.3
    x = barplot(B$b, space=space, xlim=xlim , ylim=ylim, 
                border=NA, names=receiver, las=3, ylab="Posterior samples",
                cex.axis=cex, cex.names=cex, cex.lab=cex)
    y = (x[1]/2.2)*minuslog10FDR/max(minuslog10FDR)
    tmp = sapply(1:length(x), function(i) 
        draw.circle(x[i,1], B$b[i], y[i], col=col, border=col) )
    
    # legend
    legend = seq(max(y),0,length=5)
    y2x = (ylim[2]-ylim[1])/(xlim[2]-xlim[1])
    xpos = (length(B$b)+1)*(space+1)
    ypos = ylim[2] - 3*c(1:length(legend))*max(y)*y2x
    text(xpos,ylim[2]-y2x,"FDR")
    tmp = sapply(1:length(legend), function(i) {
        draw.circle(xpos, ypos[i], legend[i], col=col, border=col) 
        text(xpos+2*max(y), ypos[i], signif(10^(-legend[i]),1) )
    } )

}




