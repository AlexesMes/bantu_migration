##Functions to calculate the accuracy and precision of models

#Calculate credible interval
credible_interval <- function(x, ci_width = 0.95)
{
  ci_upper = 0.5 + (ci_width)/2
  ci_lower = 0.5 - (ci_width)/2
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(ci_lower, ci_upper))
  return(qta)
}

#Accuracy 
accuracy <- function(sim_values, credible_interval_df){
  in_ci_range <- (sim_values >= credible_interval_df[1,] & sim_values <= credible_interval_df[2,])
  
  return(sum(in_ci_range)/length(in_ci_range))
}

#Precision
precision <- function(sim_values, credible_interval_df){
  #sim_date_range <- max(sim_values)-min(sim_values)
  average_ci_width <- mean(credible_interval_df[2,]-credible_interval_df[1,])
  
  return(average_ci_width) #sim_date_range
}
