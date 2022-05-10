#' plot_interaction_proportion
#'
#' @param prop_data a data.frame representing the proportion of types of cell-cell interaction
#' @param coordinates a data.frame representing the cell coordinates with at least two columns of "pixel_x" and "pixel_y"
#' @param x_scale a numeric value representing the scale of image zooming
#' @param cell_types_all a character vector representing the cell types to plot
#' @param img_path a character value representing the image path
#' @param scatterpie_alpha a numeric value representing the pie transparency
#' @param pie_scale a numeric value representing the pie size
#'
#' @return
#' @export
#'
#' @examples
#' 
#' setwd(system.file(package = "spacia"))
#' img_path = "data/MouseLiverST_C1/CN73_Liver_HE_C1_0.1.jpg"
#' coordinates = read.delim("data/MouseLiverST_C1/spot_ST_CN73_Liver_C1_0.1.tsv.gz", sep="\t", header=T)
#' rownames(coordinates) = with(coordinates, paste(x,y,sep="x"))
#' 
#' #' # estimate the spot radius
#' spot_radius <- calculate_spot_radius(coordinates, fct=0.25)
#' 
#' # simulate the interactive cells using the adjacent spots
#' dist <- calculate_cell_dist(coordinates, max_dist=100*spot_radius)
#' dist = subset(dist, j>i)
#' set.seed(123)
#' dist$type = sample( c("C1","C2","C3"),nrow(dist),replace=T)
#' interaction = with(dist, data.frame( i=i, j=j, weight=d, type=type ) )
#' 
#' # caculate the proportion of cell-cell interaction types
#' prop_data = t( apply( table( interaction$j, interaction$type ), 1, function(x) {
#'     if(sum(x)>0) {
#'       x/sum(x) 
#'     } else {
#'       rep(0,length(x))
#'     }
#'   } ))
#' rownames(prop_data) = rownames(coordinates)[as.integer(rownames(prop_data))]
#' 
#' cell_types_all = colnames(prop_data)
#' 
#' # plot cell-cell interaction prop
#' plot_interaction_proportion(prop_data, coordinates, x_scale=1, cell_types_all, 
#'                       img_path, scatterpie_alpha=1, pie_scale = 0.5)
#' 
plot_interaction_proportion <- function( prop_data, coordinates, x_scale=1, cell_types_all, img_path=NULL, scatterpie_alpha=1, pie_scale = 0.8 )
{
    
    library(ggplot2)
    library(scatterpie)
    
    if(0)
    {
        slice <- names(se_obj@images)[1]
        x_scale = se_obj@images[[slice]]@scale.factors$lowres
        coordinates = data.frame(se_obj@images[[slice]]@coordinates)
        prop_data = se_obj@meta.data
    }
    
    prop_data = data.frame(barcodeID=rownames(prop_data), prop_data)
    coordinates = data.frame(barcodeID=rownames(coordinates), coordinates)
    spatial_coord = merge(coordinates, prop_data)
    spatial_coord$pixel_y_scaled = spatial_coord$pixel_x*x_scale
    spatial_coord$pixel_x_scaled = spatial_coord$pixel_y*x_scale
    
    if( !is.null(img_path) ) {
        
        
        img_frmt <- base::tolower(stringr::str_sub(img_path, -4, -1))
        
        if (img_frmt %in% c(".jpg", "jpeg")) {
            img <- jpeg::readJPEG(img_path)
        } else if (img_frmt == ".png") {
            img <- png::readPNG(img_path)
        }
        
        img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
        scatterpie_plt <- suppressMessages(ggplot2::ggplot() + 
                                               ggplot2::annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) + 
                                               scatterpie::geom_scatterpie(data = spatial_coord, 
                                                                           ggplot2::aes(x = pixel_y_scaled, y = pixel_x_scaled), 
                                                                           cols = cell_types_all, color = NA, 
                                                                           alpha = scatterpie_alpha, pie_scale = pie_scale) + 
                                               ggplot2::scale_y_reverse() + 
                                               ggplot2::ylim(nrow(img), 0) + 
                                               ggplot2::xlim(0, ncol(img)) + 
                                               cowplot::theme_half_open(11, rel_small = 1) + 
                                               ggplot2::theme_void() + 
                                               ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on"))
        scatterpie_plt
        
        
            
    } else {
        
        spatial_coord$pixel_x_adjusted = max(spatial_coord$pixel_x) + min(spatial_coord$pixel_x) - spatial_coord$pixel_x
        spatial_coord$pixel_y_adjusted = max(spatial_coord$pixel_y) + min(spatial_coord$pixel_y) - spatial_coord$pixel_y
        
        scatterpie_plt <- suppressMessages(ggplot2::ggplot() + 
                                               scatterpie::geom_scatterpie(data = spatial_coord, 
                                                                           ggplot2::aes(x = pixel_x_adjusted, y = pixel_y_adjusted), 
                                                                           cols = cell_types_all, color = NA, # if NA, then colors set by default
                                                                           alpha = scatterpie_alpha, pie_scale = pie_scale) + 
                                               ggplot2::scale_y_reverse() + 
                                               ggplot2::xlim( range(spatial_coord$pixel_x_adjusted)*c(0.99, 1.01) ) + 
                                               ggplot2::ylim( range(spatial_coord$pixel_y_adjusted)*c(0.99, 1.01) ) + 
                                               cowplot::theme_half_open(11, rel_small = 1) + 
                                               ggplot2::theme_void() + 
                                               ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on"))
        
        
    }
    
    scatterpie_plt
    
}



