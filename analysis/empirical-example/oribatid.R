#### Load required libraries and helper functions
library(here)
library(ggplot2)
library(data.table)
library(vegan)
library(ade4)


source(here("R", "diffusion_map.R"))
source(here("R", "horn_overlap.R"))

#### Format data----------------------------------------------------------------
data(oribatid)
oribatid

cm <- oribatid$fau

# Convert to proportional abundance
cm <- decostand(cm, method = "total")

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

# Flip axes to for easier visualization
# site.dim[, dim1 := dim1*-1]
# site.dim[, dim2 := dim2*-1]

site.dim <- cbind(site.dim, oribatid$envir, oribatid$xy)

# Plot diffmap weighted by eigenvalues
ggplot(site.dim, aes(x = dim1, y = dim2)) + 
  geom_text(aes(label = site)) + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  theme_bw() 

ggplot(site.dim, aes(x = dim1, y = dim2, color = substrate)) +
  geom_text(aes(label = site)) + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  theme_bw() 

ggplot(site.dim, aes(x = dim1, y = dim2, color = water)) +
  geom_text(aes(label = site)) + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  theme_bw() 

#### NMDS-----------------------------------------------------------------------
hm <- HornMat(cm)
d <- 1 - hm
nm <- metaMDS(d, k = 2)
plot(nm, type = "t")

ld <- d[lower.tri(d, diag = F)]
horn.dist <- data.table(cbind(ld, dist(site.dim[, .(x, y)])))
setnames(horn.dist, c("ld", "V2"), c("horn.dist", "geo.dist"))

ggplot(horn.dist, aes(x = geo.dist, y = horn.dist)) + 
  geom_point() +
  theme_bw()


#### Diffusion distance and environmental analysis------------------------------
diff.dist <- data.table(cbind(dist(site.dim[, .(dim1, dim2)]), dist(site.dim[, .(x, y)])))
setnames(diff.dist, c("V1", "V2"), c("diff.dist", "geo.dist"))

ggplot(diff.dist, aes(x = geo.dist, y = diff.dist)) + 
  geom_point() +
  theme_bw()
