Lap <- function(sim, cmax){
  
  # This function takes a similarity matrix sim and thresholds it keeping
  # cmax nearest neighbors.
  #
  # Takes the following input:
  #
  # sim: a symmetric matrix of pairwise similarities for all samples
  #
  # cmax: an integer specifying the number of nearest neighbors to keep
  
  ## Get adjacency matrix with specified number of connections 
  adj <- sim
  
  # Set all entries above threshold to 0
  for (i in 1:nrow(adj)){
    adj[i, rank(-adj[i, ], ties.method = "min") > cmax] <- 0 
  }
  # Use ties method "min". See entry 054 for notes.
  
  # Make matrix symmetric by filling in missing values
  for(i in 1:ncol(adj)){
    for(j in 1:nrow(adj)){
      if(adj[i,j] != adj[j, i]){
        adj[i,j] <- max(adj[i,j], adj[j,i])
      }
    }
  }
  
  ## Calculate row-normalized Laplacian
  l <- -adj
  rs <- rowSums(adj)
  l <- l/rs
  diag(l) <- 1
  
  ## Return matrices as list
  out <- list(adj, l)
  names(out) <- c("adjacency", "laplacian")
  return(out)
  
}