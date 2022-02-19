#' plot_interaction_direction
#'

#' @param img_path a character value indicating the image path
#' @param coordinates a data.frame for cell coordinates, where each row represents a cell, with pixel_x and pixel_y indicating the coordinates
#' @param interaction a data frame for cell-cell interaction
#' @param x_scale a numeric value indicating the augmentation scale value relative to original image
#' @param sender_col a character value indicating the color of sender cells
#' @param reciever_col a character value indicating the color of reciever cells
#' @param interaction_col a character value indicating the color of arrows
#' @param image_alpha a numeric value indicating the background image transparency
#' @param cell_alpha a numeric value indicating the size of cells
#' @param interaction_alpha a numeric value indicating the width of interaction arrows
#' @param interaction_size a numeric value indicating the length of interaction arrows
#'
#' @return
#' @export
#'
#' @examples
#' 
#' setwd(system.file(package = "SCIA-MIL"))
#' img_path = "data/MouseLiverST_C1/CN73_Liver_HE_C1_0.1.jpg"
#' coordinates = read.delim("data/MouseLiverST_C1/spot_ST_CN73_Liver_C1_0.1.tsv.gz", sep="\t", header=T)
#' rownames(coordinates) = with(coordinates, paste(x,y,sep="x"))
#' 
#' dist <- calculate_cell_dist(coordinates, 55)
#' set.seed(123)
#' dist = dist[sample(1:nrow(dist),500), ]
#' dist = subset(dist, !j%in%i)
#' interaction = with(dist, data.frame( i=i, j=j, weight=runif(nrow(dist)) ) )
#' 
#' # plot cell-cell interaction
#' plot_interaction_direction(img_path, interaction, coordinates,x_scale=1, 
#'                            sender_col='darkred', reciever_col='steelblue', interaction_col='darkgreen',
#'                            image_alpha=0.5, cell_alpha=1, interaction_alpha=5, interaction_size=0.1)
#' 
#' 
#' 
plot_interaction_direction <- function( img_path, interaction, coordinates, x_scale=1, 
                                        sender_col='darkred', reciever_col='steelblue', interaction_col='darkgreen',
                                        image_alpha=0.5, cell_alpha=1, interaction_alpha=5, interaction_size=0.1)
{
    library(EBImage)
    
    if(image_alpha>1) image_alpha = 1
    if(image_alpha<0) image_alpha = 0
    if(cell_alpha<0) cell_alpha=0
    if(interaction_alpha<0) interaction_alpha=0
    if(interaction_size<0) interaction_size=0
    
    im <- readImage(img_path)
    im2 <- (1-image_alpha)*im + image_alpha*(im*0+1)
    
    display(im2, method='raster')
    
    weight = with(interaction, (weight-min(weight)+0.1)/(max(weight)-min(weight)+0.1) )
    
    for( k in 1:nrow(interaction) )
    {
        i = interaction[k,1]
        j = interaction[k,2]
        
        xoff=0
        yoff=0 
        xi = ( as.numeric(coordinates[i,"pixel_x"]) - xoff )*x_scale
        yi = ( as.numeric(coordinates[i,"pixel_y"]) - yoff )*x_scale
        xj = ( as.numeric(coordinates[j,"pixel_x"]) - xoff )*x_scale
        yj = ( as.numeric(coordinates[j,"pixel_y"]) - yoff )*x_scale
        
        wij = weight[k]
        #cij = interaction$type[k]
        
        points( xj, yj, pch=16, col=reciever_col, cex=cell_alpha )
        arrows(xi , yi , xj, yj, col=interaction_col, length=interaction_size, lwd=interaction_alpha*wij )
        points( xi, yi, pch=16, col=sender_col, cex=cell_alpha )
        
    }
    
}


calculate_cell_dist <- function(cells_on_spot, max_dist=Inf)
{
    do.call(rbind, lapply( 1:nrow(cells_on_spot), function(i) {   
        cat(i,"\n")
        d = sqrt( ( cells_on_spot[i, "pixel_x"] - cells_on_spot[, "pixel_x"])^2 + 
                      ( cells_on_spot[i, "pixel_y"] - cells_on_spot[, "pixel_y"])^2 ) 
        j = which(d<=max_dist)
        j = setdiff(j,i)
        res = NULL
        if(length(j)>0) res = data.frame(i=i,j=j,d=d[j])
        res
    }))
}

