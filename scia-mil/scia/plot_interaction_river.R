#' plot_interaction_river
#'
#' @param data a data.frame with at least three columns, source, target and weight
#' @param Aname a character value representing the source column name
#' @param Bname a character value representing the target column name
#' @param Wname a character value representing the weight column name
#' @param Acolor a vector of character values representing the colors of source
#' @param Bcolor a vector of character values representing the colors of target
#' @param title 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' data <- as.data.frame(UCBAdmissions)
#' plot_interaction_river(data, Aname="Admit", Bname="Dept", Wname ="Freq")
#' plot_interaction_river(data, Bname="Admit", Aname="Dept", Wname ="Freq")

plot_interaction_river <- function(data, 
                       Aname=colnames(data)[1], Bname=colnames(data)[2], Wname=colnames(data)[3], 
                       Acolor = NULL, Bcolor = NULL, 
                       title = NULL ) 
{
    
    library(ggplot2)
    library(ggalluvial)
    
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
    col.direction = "backward"
    
    gg <- ggplot( data.lodes, aes(x = factor(x, levels = ABname), y = get(Wname), 
                                              stratum = stratum, alluvium = connection, fill = stratum, label = stratum)) + 
        geom_flow(width = 1/3, aes.flow = col.direction) + 
        geom_stratum(width = 1/3, size = 0.1, color = "black", alpha = 0.8, linetype = 1) + 
        geom_text(stat = "stratum", size = font.size) + 
        scale_x_discrete(limits = c(), labels = ABname) + 
        scale_fill_manual(values = alpha(ABcolor, alpha = 0.8), drop = FALSE) + 
        theme_bw() + 
        theme(legend.position = "none", axis.title = element_blank(), axis.text.y = element_blank(), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 10)) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        ggtitle(title)

    gg
    
}





