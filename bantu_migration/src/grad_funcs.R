##Functions to plot the magnitude and direction of gradients

#Extract gradient distribution information
extract_gradinfo <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^nabla\\[',colnames(tmp))]
  return(tmp2)
}

#Determine proportion of distribution which is positive
prop_gthan_zero <- function(data) {
  # Count the number of elements greater than zero
  count_positive <- sum(data > 0)
  # Calculate the proportion
  proportion <- count_positive / length(data)
  return(proportion)
}

# Define a custom function to plot arrows
plot_arrows <- function(edges, ...) {
  for (i in 1:nrow(edges)) {
    
    ##Determine which is the start node and which is the end node
    if(edges$mean_gradient[i] >= 0){
      x_start = edges$region1_x[i]
      y_start = edges$region1_y[i]
      x_end = edges$region2_x[i]
      y_end = edges$region2_y[i] } else {
        x_start = edges$region2_x[i]
        y_start = edges$region2_y[i]
        x_end = edges$region1_x[i]
        y_end = edges$region1_y[i]
      }
    
    ##Determine the angle of the edge
    #Calculate the differences in x and y coordinates
    delta_x <- x_end - x_start
    delta_y <- y_end - y_start
    
    # Calculate the angle in radians using the arctangent function (atan2)
    angle_rad <- atan2(delta_y, delta_x)
    
    # Define the length of the arrow
    arrow_length <- abs(edges$mean_gradient[i])*40 #TODO: relative magnitude of gradient is what is important, '40' is merely a scaling parameter
    
    # Calculate the coordinates of the arrow head
    arrow_head_x <- x_end - arrow_length * cos(angle_rad)
    arrow_head_y <- y_end - arrow_length * sin(angle_rad)
    
    # Define alpha transparency value (0 to 1) depending on uncertainty in gradient
    alpha <- edges$uncertainty[i] #the proportion of the distribution in this direction
    
    arrows(x0 = x_start, y0 = y_start, x1 = arrow_head_x, y1 = arrow_head_y, col = rgb(0, 0, 1, alpha), ...)
  }
}