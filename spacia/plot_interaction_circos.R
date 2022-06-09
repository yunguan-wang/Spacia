#' plot_interaction_circos
#'
#' @param interaction a data.frame with at least three columns, source, target and weight
#' @param cell.order a vector of character values to rank the cell types in the plot
#' @param color.use a vector of colors with the same length of unique cell types
#' @param title a character value indicating the title of the plot
#'
#' @return
#' @export
#'
#' @examples
#' 
#' cells = c("C1","C2","C3","C4","C5")
#' interaction <- data.frame(source=sample(cells,20,replace=T), 
#'                           target=sample(cells,20,replace=T), 
#'                           weight=runif(20))
#' plot_interaction_circos(interaction, cell.order=cells, title="Interaction Circos plot")
#' 
#' 
plot_interaction_circos <- function(interaction, cell.order=NULL, color.use=NULL, title=NULL)
{
    library(circlize)
    #library(ComplexHeatmap)
    
    stopifnot( all(c("source", "target", "weight") %in% colnames(interaction)) )
    
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
    
    interaction <- subset(interaction, weight > 0)
    
    cell.levels <- union( as.character(interaction$source), as.character(interaction$target) )
    if( !is.null(cell.order) ) cell.levels = cell.order
    
    if (is.null(color.use)) {
        color.use <- scPalette(length(cell.levels))
        names(color.use) <- cell.levels
    } 
    
    grid.col <- color.use
    names(grid.col) <- cell.levels
    edge.color <- color.use[as.character(interaction$source)]
    
    circos.clear()
    
    chordDiagram(interaction, order = cell.levels, col = edge.color, 
                 grid.col = grid.col, transparency = 0.4, link.border = NA, 
                 directional = 1, direction.type = c("diffHeight", "arrows"), 
                 link.arr.type = "big.arrow", annotationTrack = "grid", 
                 annotationTrackHeight = 0.03, 
                 preAllocateTracks = list(track.height = max(strwidth(cell.levels))), 
                 small.gap = 1, big.gap = 10, link.visible = TRUE, 
                 scale = FALSE, group = NULL, link.target.prop = TRUE, 
                 reduce = -1)
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        xplot = get.cell.meta.data("xplot")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                    niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8) }, 
        bg.border = NA)
    
    if (!is.null(title)) {
        text(-0, 1.02, title, cex = 1)
    }
    
    circos.clear()
    
    gg <- recordPlot()
    
    gg
    
}
