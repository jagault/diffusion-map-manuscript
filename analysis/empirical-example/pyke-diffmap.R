#### Load required libraries and helper functions
library(here)
library(ggplot2)
library(data.table)
library(vegan)
library(RColorBrewer)


source(here("R", "diffusion_map.R"))
source(here("R", "horn_overlap.R"))

#### Load data------------------------------------------------------------------
# Table 1 from Pyke at al. 2001 with key for matching to Condit et al. 2001
pyke <- fread(file = paste(here(), "/analysis/empirical-example/pyke-2001.csv", 
                           sep = ""), 
              header = T)

# Environmental data from Condit et al. 2001
env <- fread(file = paste(here(), 
                          "/analysis/empirical-example/condit-site-info.csv", 
                          sep = ""), 
             header = T, check.names = T)

# Species matrix from Condit et al. 2001
spec <- fread(file = paste(here(), 
                           "/analysis/empirical-example/condit-species-matrix.csv", 
                            sep = ""), 
              header = T)

#### Subset data from Condit to match Pyke--------------------------------------
# Clean up env
setnames(env, "site.no.", "site")
env[, V8 := NULL]

# Clean up pyke
pyke[, c("plot.name", "utm.x", "utm.y") := NULL]
setnames(pyke, "key", "site")

# Select sites in Condit that match sites in Pyke
setkey(pyke, site)
setkey(env, site)

env <- merge(env, pyke)

### Make community matrix
cm <- as.matrix(spec[, -1])
rownames(cm) <- spec[, species]
cm <- t(cm) # transpose so species = columns

# Format rownames of community matrix to match site names in env
rownames(cm) <- gsub("BCI", "B", rownames(cm))
rownames(cm) <- gsub("P", "p", rownames(cm))
rownames(cm) <- gsub("p0", "p", rownames(cm))

# Select rows of community matrix present in env
cm <- cm[rownames(cm) %in% env[, site], ]


#### Calculate community dissimilarity------------------------------------------
# Calculate the sorensen dissimilarity based on presence/absence
d <- vegdist(cm, method = "bray", binary = T)


#### Diffusion map--------------------------------------------------------------
# Convert distance matrix to similarity matrix
s <- as.matrix(d)
s <- 1 - s 
diag(s) <- 0

# Get site laplacian and eigenvectors/values
site.lap <- Lap(sim = s, cmax = 10)
site.eig <- eigen(site.lap$laplacian)

# Weight eigenvectors by their respective eigenvalues
for (i in 1:ncol(site.eig$vectors)){
  site.eig$vectors[, i] <- site.eig$vectors[, i]/site.eig$values[i]
}

# Select number of dimensions to keep
k <- 1:2
# Create vector of column names
cnames <- paste("dim", rev(k), sep = "")

# Select eigenvectors corresponding to the first k dimensions
site.dim <- data.table(site.eig$vectors[, rank(site.eig$values, 
                                               ties.method = "first") %in% 
                                          (k + 1)])

# Set column names and column order
setnames(site.dim, cnames)
setcolorder(site.dim, rev(cnames))

# Add site column
site.dim[, site := rownames(cm)]

# Add env data to site.dim
setkey(site.dim, site)
site.dim <- merge(site.dim, env)

# Flip 2nd axis for plotting
site.dim[, dim2 := dim2*-1]

### Plot eigenvecots
# Make data.table of ranked eigenvalues
eig.v <- data.table(val = site.eig$values)
setorder(eig.v, val)
eig.v <- eig.v[-1,]
eig.v[, rank := rank(val)]

# Plot
ggplot(eig.v, aes(x = rank, y = 1/val)) + 
  geom_point() + 
  labs(x = "Rank", y = expression(1/lambda)) +
  theme_bw() +
  theme(text = element_text(size = 8)) 



### Plot diffusion map
# Colors for locality
lcols <- brewer.pal(7, "Dark2")
names(lcols) <- c("BCI", "Coccoli", "Ft. Sherman", "Gamboa", "Outer Watershed",
                  "Pipeline", "Santa Rita")

# Plot colored by locality
ggplot(site.dim, aes(x = dim1, y = dim2, color = locality)) +
  geom_point(size = 3) + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  scale_color_manual(values = lcols,
                    labels = c("BCI", "Coccoli", "Ft. Sherman", "Gamboa", 
                               "Outer Watershed", "Pipeline", "Santa Rita"),
                    name = "Locality") +
  theme_bw() +
  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)

# Plot sites by coordinates and colored by locality
ggplot(env, aes(x = EW.coord, y = NS.coord, color = locality)) +
  geom_point(size = 3) + 
  scale_color_manual(values = lcols,
                     labels = c("BCI", "Coccoli", "Ft. Sherman", "Gamboa", 
                                "Outer Watershed", "Pipeline", "Santa Rita"),
                     name = "Locality") +
  theme_bw()

# Plot labeled by site
ggplot(site.dim, aes(x = dim1, y = dim2)) +
  geom_text(aes(label = site), size = 1.6) + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  theme_bw() 

#### NMDS-----------------------------------------------------------------------
nm <- metaMDS(d, k = 2)
plot(nm, type = "t")

nm.dat <- data.table(nm$points)
nm.dat[, site := rownames(nm$points)]

setkey(nm.dat, site)
nm.dat <- merge(nm.dat, env)

ggplot(nm.dat, aes(x = MDS1, y = MDS2, color = locality)) +
  geom_point(size = 3) +
  scale_color_manual(values = lcols,
                     labels = c("BCI", "Coccoli", "Ft. Sherman", "Gamboa", 
                                "Outer Watershed", "Pipeline", "Santa Rita"),
                     name = "Locality") +
  theme_bw()


#### Environmental analysis within tropical moist forest life zone--------------
# Remove sites p33-p39 which are in a different floristic zone
drop <- paste(rep("p", 9), c(31:39), sep = "")
# drop <- paste(rep("p", 7), c(33:39), sep = "")
which(rownames(cm) %in% drop)
cm <- cm[-which(rownames(cm) %in% drop), ]

# Calculate the sorensen dissimilarity based on presence/absence
d <- vegdist(cm, method = "bray", binary = T)


# Convert distance matrix to similarity matrix
s <- as.matrix(d)
s <- 1 - s 
diag(s) <- 0

# Get site laplacian and eigenvectors/values
site.lap <- Lap(sim = s, cmax = 10)
site.eig <- eigen(site.lap$laplacian)

# Weight eigenvectors by their respective eigenvalues
for (i in 1:ncol(site.eig$vectors)){
  site.eig$vectors[, i] <- site.eig$vectors[, i]/site.eig$values[i]
}

# Select number of dimensions to keep
k <- 1:2
# Create vector of column names
cnames <- paste("dim", rev(k), sep = "")

# Select eigenvectors corresponding to the first k dimensions
site.dim <- data.table(site.eig$vectors[, rank(site.eig$values, 
                                               ties.method = "first") %in% 
                                          (k + 1)])

# Set column names and column order
setnames(site.dim, cnames)
setcolorder(site.dim, rev(cnames))

# Add site column
site.dim[, site := rownames(cm)]

# Add env data to site.dim
setkey(site.dim, site)
site.dim <- merge(site.dim, env)

# Flip 2nd axis for plotting
site.dim[, dim2 := dim2*-1]

### Plot eigenvecots
# Make data.table of ranked eigenvalues
eig.v <- data.table(val = site.eig$values)
setorder(eig.v, val)
eig.v <- eig.v[-1,]
eig.v[, rank := rank(val)]

# Plot
ggplot(eig.v, aes(x = rank, y = 1/val)) + 
  geom_point() + 
  labs(x = "Rank", y = expression(1/lambda)) +
  theme_bw() +
  theme(text = element_text(size = 8)) 



### Plot diffusion map
# Colors for locality
lcols <- brewer.pal(7, "Dark2")
names(lcols) <- c("BCI", "Coccoli", "Ft. Sherman", "Gamboa", "Outer Watershed",
                  "Pipeline", "Santa Rita")

# Plot colored by locality
ggplot(site.dim, aes(x = dim1, y = dim2, color = locality)) +
  geom_point(size = 3) + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  scale_color_manual(values = lcols,
                     labels = c("BCI", "Coccoli", "Ft. Sherman", "Gamboa", 
                                "Outer Watershed", "Pipeline", "Santa Rita"),
                     name = "Locality") +
  theme_bw() + 
  xlim(-3, 3) + ylim(-3, 3)

ggplot(site.dim, aes(x = dim1, y = precip)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm")

ggplot(site.dim, aes(x = dim1, y = precip)) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  geom_smooth(method = "lm")

test1 <- lm(dim1 ~ precip, site.dim)
summary(test1)

#### NMDS-----------------------------------------------------------------------
nm <- metaMDS(d, k = 2, pc = T)
plot(nm, type = "t")

nm.dat <- data.table(nm$points)
nm.dat[, site := rownames(nm$points)]

setkey(nm.dat, site)
nm.dat <- merge(nm.dat, env)

ggplot(nm.dat, aes(x = MDS1, y = MDS2, color = locality)) +
  geom_point(size = 3) +
  scale_color_manual(values = lcols,
                     labels = c("BCI", "Coccoli", "Ft. Sherman", "Gamboa", 
                                "Outer Watershed", "Pipeline", "Santa Rita"),
                     name = "Locality") +
  theme_bw()
  

ggplot(nm.dat, aes(x = -MDS1, y = precip)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm")

ggplot(nm.dat, aes(x = -MDS1, y = precip)) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  geom_smooth(method = "lm")

test2 <- lm(-MDS1 ~ precip, nm.dat)
summary(test2)



### Distance----------
diff.dist <- dist(site.dim[, .(dim1, dim2)])
p.dist <- dist(site.dim[, .(precip)])

mantel(xdis = p.dist, ydis = diff.dist)


precip <- data.table(cbind(diff.dist, p.dist))
ggplot(precip, aes(x = p.dist, y = diff.dist)) +
  geom_point() +
  theme_bw()


p.dist2 <- dist(site.dim[rownames(cm), precip])
mantel(xdis = p.dist2, ydis = d)

precip2 <- data.table(cbind(d, p.dist2))
ggplot(precip2, aes(x = p.dist2, y = d)) +
  geom_point() + 
  theme_bw()

# geo dist
diff.dist <- dist(site.dim[, .(dim1, dim2)])
geo.dist <- dist(site.dim[, .(EW.coord, NS.coord)])
geo <- data.table(cbind(diff.dist, geo.dist))
ggplot(geo, aes(x = geo.dist, y = diff.dist)) + 
  geom_point() +
  theme_bw()

geo.dist2 <- dist(site.dim[rownames(cm), .(EW.coord, NS.coord)])
geo2 <- data.table(cbind(d, geo.dist2))
geo2
ggplot(geo2, aes(x = geo.dist2, y = d)) + 
  geom_point() +
  theme_bw()
