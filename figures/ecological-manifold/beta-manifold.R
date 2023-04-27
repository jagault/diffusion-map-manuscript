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

### Response curve function-----------------------------------------------------
betaAbund <- function(pos, mode, m.abund, 
                      alpha, gamma, range){
  
  # Define variables
  m <- mode
  A <- m.abund
  a <- alpha
  g <- gamma
  x <- pos
  r <- range
  
  b <- a/(a + g)
  d <- (b^a)*(1-b)^g
  abund.x <- ((A/d)*((x - m)/r + b)^a)*(1 - ((x-m)/r + b))^g
  return(abund.x)
}

### Generate response curves----------------------------------------------------
env <- seq(0, 1, by = 0.01)

# Curve 1
c1 <- c()
for (i in 1:length(env)){
  c1[i] <- betaAbund(pos = env[i], alpha = 4, gamma = 4, m.abund = .90,
                       mode = .40, range = 1.05)
}
c1[is.na(c1)] <- 0

# Curve 2
c2 <- c()
for (i in 1:length(env)){
  c2[i] <- betaAbund(pos = env[i], alpha = 4, gamma = 4, m.abund = .45,
                     mode = .50, range = 1.05)
}
c2[is.na(c2)] <- 0

# Curve 3
c3 <- c()
for (i in 1:length(env)){
  c3[i] <- betaAbund(pos = env[i], alpha = 4, gamma = 4, m.abund = .90,
                     mode = .60, range = 1.05)
}
c3[is.na(c3)] <- 0

#### Plot-----------------------------------------------------------------------
niches <- matrix(nrow = length(env), ncol = 3)
niches[, 1:3] <- cbind(c1, c2, c3)


# Color palette for species
color <- brewer.pal(10, "Paired")

# Plot
matplot(as.data.frame(env), as.data.frame(niches), type='l', 
        xlab='Environmental Variable', ylab='Expected Species Density',
        col = color, lwd = 3)

### Sample----------------------------------------------------------------------
abundance <-  matrix(nrow = length(env), ncol = 3)
area <- 100
for (col in 1:3){
  abundance[,col] <- rpois(101, niches[,col]*area)
}

relabundance <- t(apply(abundance, 1, function(x) x/sum(x)))

#### Plot assemblage compositions-----------------------------------------------
abundance <- niches*100
relabundance <- t(apply(abundance, 1, function(x) x/sum(x)))
# Format data
cm <- relabundance
rownames(cm) <- paste("S", 1:nrow(cm), sep = "")
colnames(cm) <- LETTERS[1:ncol(cm)]

# Set as datatable for plotting
cm.dt <- data.table(cm, keep.rownames = T)
cm.dt[, env := env]
setnames(cm.dt, "rn", "site")

# Plot
ggplot(cm.dt, aes(x = A, y = B, z = C, color = env)) + 
  theme_void() + 
  scale_color_continuous(type = "viridis", name = "Env") + 
  axes_3D() + 
  stat_3D() + 
  labs_3D(labs = c("Species A", "Species B", " Species C"),
          hjust=c(0.4,0.8,1), vjust=c(2, 2.1, -0.2), angle=c(0, 0, 90))

#### Diffusion Map--------------------------------------------------------------
# Create matrix of relative abundances
cm <- as.matrix(cm.dt[, .(A, B, C)])
rownames(cm) <- cm.dt[, site]
colnames(cm) <- c("A", "B", "C")

# Calculate distances and convert to symmetric similarity matrix
s <- HornMat(cm)
diag(s) <- 0

# Get site laplacian and perform eigen decomposition
site.lap <- Lap(s, 2)
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
site.dim[, env := env]
# site.dim[, richness := apply(cm > 0, sum, MARGIN = 1)]
setcolorder(site.dim, c("site", "dim1", "dim2", "dim3", "env"))

# Plot diffmap weighted by eigenvalues
# png(filename = paste(here(), "/figures/manifold/",
#                      "diffmap.png", sep = ""),
#     height = 4.5, width = 5.5, units = "in", res = 600)
# ggplot(site.dim, aes(x = dim1, y = dim2)) + 
#   geom_point(aes(color = env), size = 5) + 
#   labs(title = "Diffusion Map",
#        x = "Dimension 1", y = "Dimension 2") + 
#   scale_color_continuous(type = "viridis", name = "Env") + 
#   theme_bw() + 
#   theme(legend.key.size = unit(0.3, "in"),
#         text = element_text(size = 30)) + 
#   xlim(-35, 35) + ylim(-35, 35)
# dev.off()

ggplot(site.dim, aes(x = dim1, y = dim2)) + 
  geom_point(aes(color = env), size = 5) + 
  labs(title = "Diffusion Map",
       x = "Dimension 1", y = "Dimension 2") + 
  scale_color_continuous(type = "viridis", name = "Env") + 
  theme_bw() + 
  theme(legend.key.size = unit(0.3, "in"),
        text = element_text(size = 30)) + 
  xlim(-350, 350) + ylim(-350, 350)

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

# Horn distance
l.hd <- hd[lower.tri(hd, diag = F)]
hd.env <- data.table(cbind(l.hd, env.dist))

ggplot(hd.env, aes(x = env.dist, y = nmds.dist)) + 
  geom_point() + 
  theme_bw() +
  labs(x = "Environmental Distance", y = "Horn Distance")