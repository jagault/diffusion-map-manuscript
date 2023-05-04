#### Load required libraries and helper functions
library(here)
library(ggplot2)
library(data.table)
library(vegan)
library(ade4)
library(viridis)


source(here("R", "diffusion_map.R"))
source(here("R", "horn_overlap.R"))

#### Format data----------------------------------------------------------------
data(doubs)
doubs

# Species matrix
cm <- doubs$fish

# Drop empty row 8
cm <- cm[-8, ]

# Convert to proportional abundance
cm <- decostand(cm, method = "total")

# Divide distance from source by 10 (see help file for Doubs)
env <- data.table(doubs$env[-8, ])
env <- cbind(env, doubs$xy[-8, ])
env[, dfs := dfs/10]

### PCA of environmnetal data---------------------------------------------------
env.pc <- rda(env, scale = T)
summary(env.pc)

biplot(env.pc, scaling = "sites")
# Loadings and biplot show that as distance from source increases, the river
# goes from oxygen rich and nutrient poor to oxygen poor and nutrient rich.

env <- cbind(env, scores(env.pc)$sites)


### Plot river course-----------------------------------------------------------
png(filename = paste(here(), "/figures/empirical-example/",
                     "river.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(env, aes(x = x, y = y)) +
  geom_path(linetype = 3) +
  geom_point(aes(color = dfs), size = 3) +
  scale_color_viridis(option = "plasma", name = "DFS") +
  labs(x = "EW Coordinate (km)", y = "NS Coordinate (km)") +
  theme_bw() +
  theme(legend.key.size = unit(0.2, "in"),
        text = element_text(size = 15))
dev.off()

#### Diffusion map--------------------------------------------------------------
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

# site.dim[, richness := apply(cm > 0, sum, MARGIN = 1)]
setcolorder(site.dim, c("site", "dim1", "dim2", "dim3"))
site.dim <- cbind(site.dim, env)

### Plot eigenvecots
# Make data.table of ranked eigenvalues
eig.v <- data.table(val = site.eig$values)
setorder(eig.v, val)
eig.v <- eig.v[-1,]
eig.v[, rank := rank(val)]

# Plot
png(filename = paste(here(), "/figures/empirical-example/",
                     "eigenvalue-plot.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(eig.v, aes(x = rank, y = 1/val)) + 
  geom_point(size = 3) + 
  labs(x = "Rank", y = expression(1/lambda)) +
  theme_bw() +
  theme(text = element_text(size = 15)) 
dev.off()

# Plot diffmap weighted by eigenvalues
png(filename = paste(here(), "/figures/empirical-example/",
                     "diffusion-map.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(site.dim, aes(x = dim1, y = dim2, color = dfs)) + 
  geom_point(size = 3) + 
  labs(x = "Dimension 1", y = "Dimension 2") +
  scale_color_viridis(option = "plasma", name = "DFS") +
  theme_bw() +
  xlim(-3.2, 3.2) + ylim(-3.2, 3.2) +
  theme(legend.key.size = unit(0.2, "in"),
        text = element_text(size = 15))
dev.off()


#### NMDS-----------------------------------------------------------------------
hm <- HornMat(cm)
d <- 1 - hm
nm <- metaMDS(d, k = 2)

site.dim <- cbind(site.dim, nm$points)

png(filename = paste(here(), "/figures/empirical-example/",
                     "nmds.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(site.dim, aes(x = MDS1, y = MDS2, color = dfs)) +
  geom_point(size = 3) +
  scale_color_viridis(option = "plasma", name = "DFS") +
  labs(x = "MDS1", y = "MDS2") +
  theme_bw() +
    theme(legend.key.size = unit(0.2, "in"),
          text = element_text(size = 15))
dev.off()

#### Distance-------------------------------------------------------------------
diff.dist <- dist(site.dim[, .(dim1)])
env.dist <- dist(site.dim[, .(dfs)])
horn.dist <- d[lower.tri(d, diag = F)]

dist.dat <- data.table(cbind(diff.dist, horn.dist, env.dist))

png(filename = paste(here(), "/figures/empirical-example/",
                     "diffusion-distance-plot.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(dist.dat, aes(x = env.dist, y = diff.dist)) +
  geom_point() +
  theme_bw() +
  labs(x = "Distance From Source (km)", y = "Diffusion Distance") +
  theme(text = element_text(size = 15))
dev.off()

png(filename = paste(here(), "/figures/empirical-example/",
                     "horn-distance-plot.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(dist.dat, aes(x = env.dist, y = horn.dist)) +
  geom_point() +
  theme_bw() +
  labs(x = "Distance From Source (km)", y = "Horn Distance") +
  theme(text = element_text(size = 15))
dev.off()