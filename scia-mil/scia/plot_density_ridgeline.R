#' plot_density_ridgeline
#'
#' @param data a data.frame of data for plotting
#' @param x_colname a character value indicating the data column name for x axis
#' @param y_colname a character value indicating the data column name for y axis
#' @param col a vector of character values indicating the colors of density plot, by default is null
#' @param title a character value indicating the title of plot
#'
#' @return
#' @export
#'
#' @examples
#' plot_density_ridgeline(data=iris, x_colname='Sepal.Length', y_colname='Species' )
#' 
#' col = c("#00AFBB", "#E7B800", "#FC4E07")
#' plot_density_ridgeline(data=iris, x_colname='Sepal.Length', y_colname='Species', col=col )
#' 
#' 
#' 
plot_density_ridgeline = function(data, 
                                  x_colname, y_colname, col=NULL, 
                                  title=paste( x_colname, "wrt.", y_colname ) )
{
    # https://www.datanovia.com/en/blog/elegant-visualization-of-density-distribution-in-r-using-ridgeline/
    
    library(ggplot2)
    library(ggridges)
    #install.packages("ggridges")
    theme_set(theme_minimal())
    
    if(is.null(col))
    {
        ggplot(data, 
               aes(x = get(x_colname), y = get(y_colname), fill = stat(x)) ) + 
            geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
            labs(x=x_colname, y=y_colname, title = title) +
            scale_fill_viridis_c(name = x_colname, option = "C")
    } else {
        
        ggplot(data, aes(x = get(x_colname), y = get(y_colname) )) +
            geom_density_ridges(aes(fill = get(y_colname))) +
            labs(x=x_colname, y=y_colname, title = title) +
            guides(fill=guide_legend(x_colname)) +
            scale_fill_manual(values = col) +
            theme(legend.position = "none")
    }
    
}


    
    