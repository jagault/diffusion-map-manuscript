# Contains a simple function for calculating pairwise differences between values
# in a vector where these values represent some measure along a corresponding
# environmental gradient. Instead of returning all pairwise differences, this
# function calculates each unique pairwise difference. 

EnvDist <- function(vec){
  # This function takes a vector of numeric values
  # Returns a matrix of unique pairwise differences corresponding to a given 
  # step along the environmental gradient. the first row gives differences for
  # all measures that are 1 env step apart, the 2nd column gives all measures
  # that are 2 steps apart, and so on. 
  
  # Create matrix to store values
  mat <- matrix(nrow = length(vec), ncol = length(vec))
  
  # Calculate distances
  for (i in 1:length(vec)){
    mat[i, ] <- abs(vec[i] - vec[i+1:length(vec)])
  }
  
  return(mat)
}

# # Short example 
# library(data.table)
# library(ggplot2)
# 
# comm <- c(5, 4, 3, 2, 1)
# env <- c(0.1, 0.2, 0.3, 0.4, 0.5)
# 
# comm.diff <- EnvDist(comm)
# env.diff <- EnvDist(env)
# 
# comm.dt <- data.table(comm.diff)
# comm.dt <- melt(comm.dt)
# comm.dt <- na.omit(comm.dt)
# setnames(comm.dt, "value", "comm.dist")
# 
# env.dt <- data.table(env.diff)
# env.dt <- melt(env.dt)
# env.dt <- na.omit(env.dt)
# setnames(env.dt, "value", "env.dist")
# 
# setkey(comm.dt, variable)
# setkey(env.dt, variable)
# dt <- cbind(comm.dt[, .(comm.dist)], env.dt[, .(env.dist)])
# 
# ggplot(dt, aes(x = env.dist, y = comm.dist)) + 
#   geom_point()