#### Load required libraries and helper functions
library(here)
library(data.table)
library(vegan)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(gg3D)

source(here("R", "diffusion_map.R"))
source(here("R", "horn_overlap.R"))
### Generate response curves----------------------------------------------------
set.seed(10)

M <- 3 # size of species inventory
N <- 25 # sample size
nichepos <- seq(0.1, 0.9, length.out = M)
nichewidth <- runif(M, 0.4, 0.5)

xseq <- seq(0, 1, by = 0.01)
lambda <- function(x, pos, width){1/width*exp(-(x-pos)^2/(width^2))} # Gaussian bell curve as niche profile
niches <- matrix(nrow = length(xseq), ncol = M)
for(col in 1:M){
  niches[,col] <- lambda(xseq, nichepos[col], nichewidth[col])
}

#### Plot response curves-------------------------------------------------------
# Color palette for species
color <- brewer.pal(10, "Paired")

# Plot
# png(filename = paste(here(), "/figures/manifold/",
#                      "repsonse-curves.png", sep = ""),
#     height = 4.5, width = 5.5, units = "in", res = 600)
matplot(as.data.frame(xseq), as.data.frame(niches), type='l', 
        xlab='Environmental Variable', ylab='Expected Species Density',
        col = color, lwd = 3)
# dev.off()

#### Sample from response curves to generate assemblages------------------------
xsample <- seq(0, 1, length.out = N) # evenly spaced environmental samples
dens <- matrix(nrow = N, ncol = M)
abundance <-  matrix(nrow = N, ncol = M)
area <- 20

# Multiple realizations for each environmental value
n.sim <- 10 # Number of realizations
relabundance <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  
  for (col in 1:M){
    dens[,col] <- lambda(xsample,  nichepos[col], nichewidth[col])
    # abundance[,col] <- round(dens[,col]*area)
    abundance[,col] <- rpois(N, dens[,col]*area)
  }
  relabundance[[i]] <- data.table(t(apply(abundance, 1, function(x) x/sum(x))))
}


#### Plot assemblage compositions-----------------------------------------------
# Format data
cm.dt <- rbindlist(relabundance, idcol = T)
setnames(cm.dt, c(".id", "V1", "V2", "V3"), c("sample", "A", "B", "C"))
cm.dt[, site := rep(paste("S", 1:N, sep = ""), times = n.sim)]
cm.dt[, env := rep(xsample, times = n.sim)]

# Plot
# png(filename = paste(here(), "/figures/manifold/",
#                      "manifold.png", sep = ""),
#     height = 7, width = 7, units = "in", res = 600)
ggplot(cm.dt, aes(x = A, y = B, z = C, color = env)) + 
  theme_void() + 
  scale_color_continuous(type = "viridis", name = "Env") + 
  axes_3D() + 
  stat_3D() + 
  labs_3D(labs = c("Species A", "Species B", " Species C"),
          hjust=c(0.4,0.8,1), vjust=c(2, 2.1, -0.2), angle=c(0, 0, 90))
# dev.off()


#### Diffusion Map--------------------------------------------------------------
# Create matrix of relative abundances
cm.dt[, site.samp := paste(site, sample, sep = ".")]
cm <- as.matrix(cm.dt[, .(A, B, C)])
rownames(cm) <- cm.dt[, site.samp]
colnames(cm) <- c("A", "B", "C")

# Calculate distances and convert to symmetric similarity matrix
s <- HornMat(cm)
diag(s) <- 0

# Get site laplacian and perform eigen decomposition
site.lap <- Lap(s, 10)
site.eig <- eigen(site.lap$laplacian)

# Weight eigenvector values by eigenvalue
for (i in 1:ncol(site.eig$vectors)){
  site.eig$vectors[, i] <- site.eig$vectors[, i]/site.eig$values[i]
}

# Select number of dimensions to keep
k <- 1:3
# Create vector of column names
cnames <- paste("dim", rev(k), sep = "")

# Select eigenvectors corresponding to the first k dimensions
site.dim <- data.table(site.eig$vectors[, rank(site.eig$values, 
                                               ties.method = "first") %in% 
                                          (k + 1)])
# Set column names and column order
setnames(site.dim, cnames)
setcolorder(site.dim, rev(cnames))

# Add site and env column
site.dim[, site := rownames(cm)]
site.dim[, site:= as.factor(site)]
site.dim[, env := rep(xsample, times = n.sim)]
# site.dim[, richness := apply(cm > 0, sum, MARGIN = 1)]
setcolorder(site.dim, c("site", "dim1", "dim2", "dim3", "env"))

# Plot diffmap weighted by eigenvalues
# png(filename = paste(here(), "/figures/manifold/",
#                      "diffmap.png", sep = ""),
#     height = 4.5, width = 5.5, units = "in", res = 600)
ggplot(site.dim, aes(x = dim1, y = dim2)) + 
  geom_point(aes(color = env), size = 5) + 
  labs(title = "Diffusion Map",
       x = "Dimension 1", y = "Dimension 2") + 
  scale_color_continuous(type = "viridis", name = "Env") + 
  theme_bw() + 
  theme(legend.key.size = unit(0.3, "in"),
        text = element_text(size = 30)) + 
  xlim(-35, 35) + ylim(-35, 35)
# dev.off()

ggplot(site.dim, aes(x = dim1, y = dim2)) + 
  geom_point(aes(color = env), size = 5) + 
  labs(title = "Diffusion Map",
       x = "Dimension 1", y = "Dimension 2") + 
  scale_color_continuous(type = "viridis", name = "Env") + 
  theme_bw() + 
  theme(legend.key.size = unit(0.3, "in"),
        text = element_text(size = 30)) 

#### NMDS-----------------------------------------------------------------------
hs <- HornMat(cm)
hd <- 1 - hs
tol <- 1e-10
hd[hd <= tol] <- 0

test <- metaMDS(hd, k = 2)
plot(test, type = "t")

test$points

#### Distance-------------------------------------------------------------------
# Environmental distances
env.dist <- round(dist(site.dim[, .(env)]), 4)

# Diffusion distance
diff.dist <- round(dist(as.matrix(site.dim[, .(dim1, dim2)])), 4)
diff.env <- data.table(cbind(diff.dist, env.dist))

diff.cor <- diff.env[, .(spearman = cor(diff.dist, env.dist, 
                                          method = "spearman"),
                           pearson = cor(diff.dist, env.dist))]

ggplot(diff.env, aes(x = env.dist, y = diff.dist)) + 
  geom_point() + 
  theme_bw() +
  labs(x = "Environmental Distance", y = "Diffusion Distance")

# NMDS distance
nmds.dist <- round(dist(test$points), 4)
nmds.env <- data.table(cbind(nmds.dist, env.dist))

nmds.cor <- nmds.env[, .(spearman = cor(nmds.dist, env.dist, 
                                        method = "spearman"),
                         pearson = cor(nmds.dist, env.dist))]

ggplot(nmds.env, aes(x = env.dist, y = nmds.dist)) + 
  geom_point() + 
  theme_bw() +
  labs(x = "Environmental Distance", y = "NMDS Distance")

# Another try
test.dat <- data.table(test$points)
test.dist <-  round(dist(as.matrix(test.dat[, .(MDS1, MDS2)])), 4)
test.env <- data.table(cbind(test.dist, env.dist))

test.cor <- test.env[, .(spearman = cor(test.dist, env.dist, 
                                        method = "spearman"),
                         pearson = cor(test.dist, env.dist))]

ggplot(test.env, aes(x = env.dist, y = test.dist)) + 
  geom_point() + 
  theme_bw() +
  labs(x = "Environmental Distance", y = "NMDS Distance")

# Diff dimension 1
# Diffusion distance
diff.dist <- round(dist(as.matrix(site.dim[, .(dim1)])), 4)
diff.env <- data.table(cbind(diff.dist, env.dist))

diff.cor <- diff.env[, .(spearman = cor(diff.dist, env.dist, 
                                        method = "spearman"),
                         pearson = cor(diff.dist, env.dist))]

ggplot(diff.env, aes(x = env.dist, y = diff.dist)) + 
  geom_point() + 
  theme_bw() +
  labs(x = "Environmental Distance", y = "Diffusion Distance")
