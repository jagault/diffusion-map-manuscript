GetLap <- function(mat, cmax){
  
  # This function takes a matrix of objects as rows with associated data as 
  # columns, calculates the similarity between objects, converts the similarity
  # matrix to an adjacency matrix with number of connections between nodes 
  # (objects) specified by cmax, and calculates the corresponding Laplacian. 
  # Returns the similarity, adjacency, and Laplacian matrices as a list.
  #
  #
  # Takes the following input:
  # mat: an array with entries of class numeric. 
  #
  # cmax: a number specifying the threshold number of connections between nodes
  #
  #
  # Returns: a list with the following entries
  # $similarity: a matrix containing the similarities between objects
  #
  # $adjacency: a symmetric, weighted adjacency matrix of objects where the 
  # number of connections between nodes (objects) is determined by cmax. 
  #
  # $laplacian: the row-normalized laplacian matrix calculated from the 
  # adjacency matrix
  
  
  ## Standardize columns
  std.mat <- apply(mat, MARGIN = 2, FUN = function(x) (x - mean(x))/sd(x))
  
  ## Calculate Euclidean distance and convert to similarity
  d <- dist(std.mat, method = "euclidean") # returns class 'dist'
  # convert to similarity
  s <- 1/d 
  # convert from class 'dist' to full matrix
  s <- as.matrix(s) 
  
  ## Get adjacency matrix with specified number of connections 
  adj <- s
  
  # Set all entries above threshold to 0
  for (i in 1:nrow(adj)){
    adj[i, rank(-adj[i, ], ties.method = "random") > cmax] <- 0 
  }
  
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
  out <- list(s, adj, l)
  names(out) <- c("similarity", "adjacency", "laplacian")
  return(out)
  
}



GetZeros <- function(mat){
  
  # This function takes a matrix of objects and, based on the euclidean distance
  # between them, gets the unique pairs of objects that are identical. If no 
  # objects are identical, returns message "No identical rows in data matrix." 
  # If there are unique pairs of objects, returns the unique compartments of 
  # connections as a list. Eeach element of the list contains the names of 
  # objects that are equal to each other. 
  
  # Calculate distance
  dist.mat <- dist(mat, method = "euclidean")
  dist.mat <- as.matrix(dist.mat)
  diag(dist.mat) <- 1 
  
  if (any(dist.mat == 0) == F){
    
    message("No identical rows in data matrix.")
    
  } else {
    
    # Get index of identical values
    dist.mat[upper.tri(dist.mat)] <- 1 
    ind <- which(dist.mat == 0, arr.ind = T)
    ind.name <- ind
    
    # Select names corresponding to index of identical values
    # This gives matrix of unique pairs
    ind.name[, 1] <- rownames(dist.mat)[ind[ ,1]]
    ind.name[, 2] <- colnames(dist.mat)[ind[, 2]]
    rownames(ind.name) <- colnames(ind.name) <- NULL
    
    # Get all connections from matrix of unique pairs
    # Create a list which will be length nrow(ind.name)
    # For each row in ind.name, create an element of the list which tells which
    # rows share at least one entry
    overlap <- list() 
    for (i in 1:nrow(ind.name)){
      overlap[[i]] <- apply(ind.name, MARGIN = 1, 
                            function(x, y) any(is.element(x, ind.name[i, ])))
    }
    
    # Get the unique entries of overlap. These are the unique compartments of 
    # connections and will serve as an index to select the entry names
    overlap <- unique(overlap)
    
    # Create a list to fill with the names of the overlapping entries. Each 
    # element of the list contains a unique compartment of connections and the 
    # names of the entries contained within. 
    connections <- list()
    for (i in 1:length(overlap)){
      connections[[i]] <- unique(as.vector(ind.name[overlap[[i]], ]))
    }
    
    return(connections)
    
  }
}


GetLapNS <- function(mat, cmax){
  
  # This function takes a matrix of objects as rows with associated data as 
  # columns, calculates the similarity between objects, converts the similarity
  # matrix to an adjacency matrix with number of connections between nodes 
  # (objects) specified by cmax, and calculates the corresponding Laplacian. 
  # Returns the similarity, adjacency, and Laplacian matrices as a list.
  # 
  # Unlike GetLap, this function does not standardize columns. 
  #
  # Takes the following input:
  # mat: an array with entries of class numeric. 
  #
  # cmax: a number specifying the threshold number of connections between nodes
  #
  #
  # Returns: a list with the following entries
  # $similarity: a matrix containing the similarities between objects
  #
  # $adjacency: a symmetric, weighted adjacency matrix of objects where the 
  # number of connections between nodes (objects) is determined by cmax. 
  #
  # $laplacian: the row-normalized laplacian matrix calculated from the 
  # adjacency matrix
  
  
  ## Calculate Euclidean distance and convert to similarity
  d <- dist(mat, method = "euclidean") # returns class 'dist'
  # convert to similarity
  s <- 1/d 
  # convert from class 'dist' to full matrix
  s <- as.matrix(s) 
  
  ## Get adjacency matrix with specified number of connections 
  adj <- s
  
  # Set all entries above threshold to 0
  for (i in 1:nrow(adj)){
    adj[i, rank(-adj[i, ], ties.method = "random") > cmax] <- 0 
  }
  
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
  out <- list(s, adj, l)
  names(out) <- c("similarity", "adjacency", "laplacian")
  return(out)
  
}

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
    adj[i, rank(-adj[i, ], ties.method = "random") > cmax] <- 0 
  }
  
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