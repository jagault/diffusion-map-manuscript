# This script contains functions to return the horn measure of overlap between
# all pairs of sites in a community matrix.
# The overlap is defined following "Horn (1966). Measures of Overlap in 
# Comparative Ecological Studies." 

HornIndex <- function(a, b) {
  # Takes two vectors containing proportional abundances
  # Returns the horn index of overlap between them
  ind <- (sum((a + b)*log(a + b), na.rm = T) - 
            sum(a*log(a), na.rm = T) - 
            sum(b*log(b), na.rm = T))/
    (2*log(2) - log(1) - log(1))
  tol <- 1e-10 # Tolerance 
  ind[ind <= tol] <- 0 # Set value less than tolerance to zero
  return(ind)
}

HornMat <- function(mat){
  # Takes a matrix where rows are samples and columns are species
  # Rows contain proportional abundance
  # Returns a matrix of all pairwise Horn overlaps
  dm <- matrix(nrow = nrow(mat),
               ncol = nrow(mat))
  rownames(dm) <- colnames(dm) <- rownames(mat)
  for (i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      dm[i,j] <- HornIndex(mat[i,], mat[j,])
    }
    
  }
  return(dm)
}


# # Full proportional overlap
# s1 <- c(10, 8, 6, 4)
# s2 <- c(10, 8, 6, 4)
# 
# # Full proportional overlap
# s1 <- c(10, 8, 6, 4)
# s2 <- c(5, 4, 3, 2)
# 
# # 50% overlap? Gives same answer as morisita-horn
# s1 <- c(1, 1, 1, 0)
# s2 <- c(1, 0, 1, 1)
# 
# # Compositional shift. Morisita-horn gives slightly lower overlap
# s1 <- c(10, 8, 6, 4)
# s2 <- c(4, 6, 8, 10)
# 
# # Another compositional shift. morisita horn-gives even lower similarity
# # The weighted, raw abundance horn index is now slightly different. This
# # is most likely due to difference in total # of individuals
# s1 <- c(10, 8, 8, 8)
# s2 <- c(10, 1, 1, 1)
# 
# # Large difference in total abundance. Both weighted and unweighted versions 
# # give the same overlap. Why is this? I'm guessing it's because they are 
# # proportionally identical which will give score of 1 no matter what. 
# # I think weighting only makes a difference if the sites are not proportionally
# # identical 
# s1 <- c(100, 80, 60, 40)
# s2 <- c(10, 8, 6, 4)
# 
# 
# # Community matrix
# cm <- rbind(s1, s2)
# colnames(cm) <- letters[1:4]
# cm
# 
# # first let's try the horn distance based on frequencies
# # transform to frequencies
# cm.p <- decostand(cm, method = "total")
# cm.p
# 
# # Horn overlap for frequencies
# a <- c()
# b <- c()
# c <- c()
# for (i in 1:ncol(cm.p)){
#   a[i] <- (cm.p[1, i] + cm.p[2, i])*log(cm.p[1, i] + cm.p[2, i])
#   b[i] <- (cm.p[1, i])*log(cm.p[1, i])
#   c[i] <- (cm.p[2, i])*log(cm.p[2, i])
# }
# 
# d <- sum(na.omit(a)) - sum(na.omit(b)) - sum(na.omit(c))
# d/(2*log(2) - log(1) - log(1))
# 
# HornMat(cm.p)
# 
# # Horn overlap for raw abundances
# a <- c()
# b <- c()
# c <- c()
# for (i in 1:ncol(cm)){
#   a[i] <- (cm[1, i] + cm[2, i])*log(cm[1, i] + cm[2, i])
#   b[i] <- (cm[1, i])*log(cm[1, i])
#   c[i] <- (cm[2, i])*log(cm[2, i])
# }
# 
# X <- sum(cm[1, ])
# Y <- sum(cm[2, ])
# 
# d <- sum(na.omit(a)) - sum(na.omit(b)) - sum(na.omit(c))
# d/((X + Y)*log(X + Y) - X*log(X) - Y*log(Y))
# 
# 1 - vegdist(cm.p, method = "horn", binary = FALSE)