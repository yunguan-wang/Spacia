#' plot_interaction_direction
#'
#' @param img_path a character value indicating the image path
#' @param coordinates a data.frame for cell coordinates, where each row represents a cell, with pixel_x and pixel_y indicating the coordinates
#' @param interaction a data frame for cell-cell interaction
#' @param xoff a numeric value indicating the offset on x-axis relative to the original image
#' @param yoff a numeric value indicating the offset on y-axis relative to the original image 
#' @param x_scale a numeric value indicating the augmentation scale value relative to original image
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
#' # estimate the spot radius
#' spot_radius <- calculate_spot_radius(coordinates, fct=0.25)
#' 
#' # simulate the interactive cells using the adjacent spots
#' dist <- calculate_cell_dist(coordinates, 4*spot_radius)
#' set.seed(123)
#' dist$type = sample(1:3,nrow(dist),replace=T)
#' dist = subset(dist, j>i)
#' interaction = with(dist, data.frame( i=i, j=j, weight=d, type=type ) )
#' 
#' # plot cell-cell interaction
#' plot_interaction_direction(img_path, interaction, coordinates)
#' 
#' 
plot_interaction_direction <- function( img_path, interaction, coordinates, xoff=0, yoff=0, x_scale=1 )
{
    library(EBImage)
    
    array_type <- "1k"
    fct <- ifelse(array_type == "1k", 0.25, 50/(sqrt(2)*100))
    
    im <- readImage(img_path)
    
    display(im,method='raster')
    plot_spot_info(coordinates, spot_cols='yellow', barcode=F, fct=fct)
    
    for( k in 1:nrow(interaction) )
    {
        i = interaction[k,1]
        j = interaction[k,2]
        
        xi = ( as.numeric(coordinates[i,"pixel_x"]) - xoff )*x_scale
        yi = ( as.numeric(coordinates[i,"pixel_y"]) - yoff )*x_scale
        xj = ( as.numeric(coordinates[j,"pixel_x"]) - xoff )*x_scale
        yj = ( as.numeric(coordinates[j,"pixel_y"]) - yoff )*x_scale
        
        wij = interaction$weight[k]
        cij = interaction$type[k]
        
        arrows(xi , yi , xj, yj, col=cij, length=0.1, lwd=wij/100 )
        
    }
    
}


calculate_spot_radius <- function(spot_coordinates, fct)
{
    distMat <- dist(spot_coordinates[, c("pixel_x", "pixel_y")])
    distMat <- as.matrix(distMat)
    diag(distMat) <- Inf
    center_to_center_distance <- mean(apply(distMat, 2, min))
    
    cat("Defining spot radius...\n")
    center_to_center_distance*fct
    
}

plot_spot_info <- function(spot_coordinates, xoff=0, yoff=0, x_scale=1, 
                           spot_cols='green', barcode=T, fct=0.25)
{
    
    if(length(spot_cols)==1) spot_cols = rep( spot_cols, nrow(spot_coordinates) )
    
    plot_circle <- function(x, y, r, c) {
        angles <- seq(0, 2*pi, length.out = 360)
        lines(r*cos(angles) + x, r*sin(angles) + y, col=c)
    }
    
    if( all(c("X","Y") %in% colnames(spot_coordinates)) )
    {
        spot_coordinates$pixel_x = (spot_coordinates$X-xoff)*x_scale
        spot_coordinates$pixel_y = (spot_coordinates$Y-yoff)*x_scale
    }
    
    if( all(c("imagecol","imagerow") %in% colnames(spot_coordinates)) )
    {
        spot_coordinates$pixel_x = (spot_coordinates$imagecol-xoff)*x_scale
        spot_coordinates$pixel_y = (spot_coordinates$imagerow-yoff)*x_scale
    }
    
    spot_radius <- calculate_spot_radius(spot_coordinates, fct)
    
    for(i in 1:nrow(spot_coordinates)) {
        
        x = as.numeric(spot_coordinates[i, "pixel_x"])
        y = as.numeric(spot_coordinates[i, "pixel_y"])
        r = spot_radius
        c = spot_cols[i]
        plot_circle(x, y, r, spot_cols)
        
        #text(x, y-r, as.character(spot_coordinates$barcode[i]) )
        if(barcode) text(x, y-r, rownames(spot_coordinates)[i])
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

