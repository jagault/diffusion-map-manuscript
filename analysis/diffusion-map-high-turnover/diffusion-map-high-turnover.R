#### Load required libraries and helper functions
library(here)
library(pracma)
library(ggplot2)
library(data.table)
library(vegan)
library(ape)

source(here("R", "diffusion_map.R"))
source(here("R", "horn_overlap.R"))

#### Generate Response Surfaces-------------------------------------------------
set.seed(450)

# Size of species inventory 
n.spec <-  100

# Initialize list to store values related to covariance matrix 
niche <- list()

# Specify mean of species response surface for env vars x and y
niche$centroid <- matrix( runif(2*n.spec, 0, 1), nrow = 2, ncol = n.spec) 

### Generate Variance/Covariance matrix
# Randomly select values for scaling factors and cos(theta)
niche$l1 <- runif(n.spec, 0.005, 0.01) # Scaling factor lambda 1
niche$l2 <- runif(n.spec, 0.005, 0.01) # Scaling factor lambda 2
niche$cos <- runif(n.spec, -1, 1) # Cos(theta)
niche$sin <- sample( c(-1, 1), replace = T, n.spec) * sqrt( 1 - niche$cos^2 ) # Sin(theta)

# Initialize array to store values for covariance matrix
niche$Cov <- array(dim = c(2,2,n.spec))

# Calculate variance and covariance values
niche$Cov[1,1,] <- niche$l1 * niche$cos^2 + niche$l2 * niche$sin^2
niche$Cov[2,2,] <- niche$l1 * niche$sin^2 + niche$l2 * niche$cos^2
niche$Cov[1,2,] <- niche$Cov[2,1,] <- (niche$l1 - niche$l2) * niche$cos * niche$sin


### Evaluate response surfaces on environmental grid
# Create environmental grid from 0-1 in steps of 0.2
x <- seq(0, 1, by = 0.02)
mesh <- meshgrid(x)
meshXY <- rbind(as.vector(mesh$X), as.vector(mesh$Y))

# Define functions for evaluating probability density function
quad.form <- function(A, x) {
  colSums(x * (A %*% x))
}

bell <- function(C, det, dist) {
  1/sqrt(det) * exp(-quad.form(C,dist))
}

# Initialize array to store values of probability density function
niche$shape <- array(dim = c(length(x),length(x),n.spec))

# Evaluate probability density function across environmental grid
# Used for contour plot below
for (spec in 1:n.spec){
  C <- solve(niche$Cov[,,spec])
  det <- niche$l1[spec]*niche$l2[spec]
  dist <- meshXY - niche$centroid[,spec]
  niche$shape[,,spec] <- matrix(bell(C, det, dist), nrow = length(x))
}


#### Contour plot of response surfaces------------------------------------------
# Add rownames and colnames
rownames(niche$shape) <- colnames(niche$shape) <- x

# Make list for converting array to data.table 
dt.list <- vector(mode = "list", length = n.spec)

# Convert each matrix of the array to a data.table
dt.list <- apply(niche$shape, 3, data.table, keep.rownames = T)

# Bind list into a data.table
dens.dt <- rbindlist(dt.list, idcol = T)
setnames(dens.dt, c(".id", "rn"), c("spec", "y"))
dens.dt <- melt(dens.dt, id.vars = c("spec", "y"), variable.name = "x", 
                value.name = "dens", variable.factor = F)

# Convert columns to correct classes
dens.dt[, spec := as.factor(spec)] 
dens.dt <- dens.dt[, lapply(.SD, as.numeric), by = spec]

# Plot
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "contour-plot.png", sep = ""),
    height = 3, width = 3, units = "in", res = 600)
ggplot(dens.dt, aes(x = x, y = y, z = dens, color = spec)) + 
  geom_contour(show.legend = F) + 
  theme_bw()
dev.off()


#### Generate community realizations--------------------------------------------
# Sample size = n.samp^2
n.samp <- 10

set.seed(20)
# Equally spaced samples along both environmental gradients.
x.samp <- seq(0.1, 0.9, length.out = n.samp)
samp.mesh <- meshgrid(x.samp)
grid.samp <- rbind(as.vector(samp.mesh$X), as.vector(samp.mesh$Y))

# Area to scale densities
area <- 10

# Initialize array to store values of probability density function
dens <- array(dim = c(1, n.samp^2, n.spec))

# Make dataframe for plotting
samp <- data.table(t(grid.samp))
samp[, site := paste("S", 1:nrow(samp), sep = "")]
setnames(samp, c("V1", "V2"), c("x.samp", "y.samp"))

# Evaluate pdf for each env grid point for each species 
comm.lambda <- matrix(nrow = n.samp^2, ncol = n.spec)
for (spec in 1:n.spec){
  C <- solve(niche$Cov[,,spec])  
  det <- niche$l1[spec]*niche$l2[spec]
  dist <- grid.samp - niche$centroid[,spec]
  dens[,,spec] <- bell(C, det, dist)
  lambda <- dens * area
  comm.lambda[, spec] <- as.vector(lambda[,,spec])
}

# Generate multiple realizations of assemblages
n.sim <- 100 # Number of realizations

comm <- vector(mode = "list", length = n.sim)
m <- matrix(nrow = n.samp^2, ncol = n.spec)
rownames(m) <- paste("site", seq(n.samp^2), sep = "")
colnames(m) <- paste("spec", seq(n.spec), sep = "")
for (i in 1:n.sim) {
  for (spec in 1:n.spec){
    m[, spec] <- rpois(n.samp^2, comm.lambda[, spec])
  }
  comm[[i]] <- m
}

#### Plot sampling grid on contour plot-----------------------------------------
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "contour-sample.png", sep = ""),
    height = 3, width = 3, units = "in", res = 600)
ggplot() +
  geom_contour(data = dens.dt, aes(x = x, y = y, z = dens, color = spec),
               show.legend = F) + 
  geom_point(data = samp, aes(x = x.samp, y = y.samp), size = 0.5) + 
  theme_bw() + 
  labs(x = "Environmental Variable 1", y = "Environmental Variable 2")
dev.off()


#### Summary of number of disjoint samples--------------------------------------
# Calculate proportional abundances
cm <- vector(mode = "list", length = n.sim)
for (i in 1:n.sim){
  cm[[i]] <- t(apply(comm[[i]], 1, function(x) x/sum(x)))
  
}


# Calculate Horn overlap
hm <- vector(mode = "list", length = n.sim)
for (i in 1:n.sim) {
  hm[[i]] <- HornMat(cm[[i]])
}

# Calculate % of sites with no shared species
dj <- vector(length = n.sim)
for (i in 1:n.sim){
  dj[i] <- sum(hm[[i]][lower.tri(hm[[i]], diag = F)] == 0)/((n.sim^2 - n.sim)/2)*100
}

dj.summary <- summary(dj)
dj.summary <- round(unclass(dj.summary), 5)
dj.summary <- data.table(t(dj.summary))

write.csv(dj.summary, file = paste(here(), 
                                   "/analysis/diffusion-map-high-turnover/",
                                   "disjoint-summary.csv", sep = ""))


#### Diffusion map and plot-----------------------------------------------------
# Calculate distances and convert to symmetric similarity matrix
s <- hm # Use horn overlap already calculated above

for (i in 1:n.sim) {
  diag(s[[i]]) <- 0
}

# Get site laplacian and perform eigen decomposition
site.lap <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  site.lap[[i]] <- Lap(s[[i]], 10)
}

site.eig <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  site.eig[[i]] <- eigen(site.lap[[i]]$laplacian)
}

# Weight eigenvectors by their respective eigenvalue
for (i in 1:n.sim) {
  for (j in 1:ncol(site.eig[[i]]$vectors)){
    site.eig[[i]]$vectors[, j] <- site.eig[[i]]$vectors[, j]/site.eig[[i]]$values[j]
  }
}

# Select number of dimensions to keep
k <- 1:2
# Create vector of column names
cnames <- paste("dim", rev(k), sep = "")

# Select eigenvectors corresponding to the first k dimensions
dim.list <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  dim.list[[i]] <- data.table(site.eig[[i]]$vectors[,rank(site.eig[[i]]$values, 
                                                          ties.method = "first") 
                                                    %in% (k + 1)])
}

# Convert to single datatable
site.dim <- rbindlist(dim.list, idcol = T)

# Set column names and column order
setnames(site.dim, c("sim", cnames))
setcolorder(site.dim, rev(cnames))

# Add site and env column
sn <- paste("site", seq(n.samp^2), sep = "")
site.dim[, site := rep(sn, length.out = nrow(site.dim))]
site.dim[, site := as.factor(site)]
site.dim[, sim := as.factor(sim)]


# Make data.table of ranked eigenvalues for diffusion map
eig.list <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  eig.list[[i]] <- data.table(val = site.eig[[i]]$values)
  setorder(eig.list[[i]], val)
  eig.list[[i]] <- eig.list[[i]][-1,]
  eig.list[[i]][, rank := rank(val)]
}

eig.v <- rbindlist(eig.list, idcol = T)
setnames(eig.v, ".id", "sim")
eig.v[, sim := as.factor(sim)]

# Eigenvalue Plot
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "eigenvalue-plot.png", sep = ""),
    height = 2, width = 3, units = "in", res = 600)
ggplot(eig.v, aes(x = rank, y = 1/val)) + 
  geom_point(size = 0.1, show.legend = F) + 
  labs(x = "Rank", y = expression(1/lambda)) +
  theme_bw()
dev.off()

### Summarize eigenvalues to plot mean and std dev
eig.sd <- copy(eig.v)
eig.sd[, val := 1/val]
eig.sd <- eig.sd[, .(mean = mean(val), sd = sd(val)), by = rank]

# Eigenvalue plot with mean and std dev
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "eigenvalue-plot-std-dev.png", sep = ""),
    height = 2, width = 3, units = "in", res = 600)
ggplot(eig.sd, aes(x = rank, y = mean)) + 
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = mean  - sd, ymax = mean + sd), width = 0.1) +
  labs(x = "Rank", y = expression(1/lambda)) +
  theme_bw()
dev.off()

# Single diffusion map
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "single-diffusion-map-plot.png", sep = ""),
    height = 3, width = 3, units = "in", res = 600)
ggplot(site.dim[sim == 1,], aes(x = dim1, y = dim2)) + 
  geom_point(show.legend = F, size = 0.5) + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  theme_bw() 
dev.off()


#### Compare diffusion and Horn distances---------------------------------------
# Add environmental variables to site.dim
site.dim[, env.x := rep(grid.samp[1, ], length.out = nrow(site.dim))]
site.dim[, env.y := rep(grid.samp[2, ], length.out = nrow(site.dim))]

env.dist <- round(dist(t(grid.samp)), 4)

dist.list <- vector(mode = "list", length = n.sim)
for (i in 1:n.sim) {
  dist.list[[i]] <- round(dist(as.matrix(site.dim[sim == i, .(dim1, dim2)])), 4)
  dist.list[[i]] <- data.table(cbind(dist.list[[i]], env.dist))
}

dist.frame <- rbindlist(dist.list, idcol = T)
setnames(dist.frame, c(".id", "V1"), c("sim", "diff.dist"))
dist.frame[, sim := as.factor(sim)]

# Calculate correlation
diff.cor <- dist.frame[, .(spearman = cor(diff.dist, env.dist, 
                                          method = "spearman"),
                           pearson = cor(diff.dist, env.dist)), by = sim]

# Summarize and write to file
# Spearman
diff.cor.spearman <- summary(diff.cor[, spearman])
diff.cor.spearman <- round(unclass(diff.cor.spearman), 5)
diff.cor.spearman <- data.table(t(diff.cor.spearman))

write.csv(diff.cor.spearman, 
          file = paste(here(), 
                       "/analysis/diffusion-map-high-turnover/",
                       "diffmap-spearman-summary.csv", 
                       sep = ""))

# Pearson
diff.cor.pearson <- summary(diff.cor[, pearson])
diff.cor.pearson <- unclass(diff.cor.pearson)
diff.cor.pearson <- data.table((t(diff.cor.pearson)))

write.csv(diff.cor.pearson, 
          file = paste(here(), 
                       "/analysis/diffusion-map-high-turnover/",
                       "diffmap-pearson-summary.csv", 
                       sep = ""))


# Plot diffusion distance vs environmental distance
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "diffusion-distance-plot.png", sep = ""),
    height = 2, width = 3, units = "in", res = 600)
ggplot(dist.frame, aes(x = env.dist, y = diff.dist)) + 
  geom_point(show.legend = F, size = 0.1) + 
  theme_bw() +
  labs(x = "Environmental Distance", y = "Diffusion Distance")
dev.off()

#### Horn Distance
hd.list <- vector(mode = "list", length = n.sim)
tol <- 1e-10 # Tolerance 

for (i in 1:n.sim) {
  hd.list[[i]] <- 1 - hm[[i]]
  hd.list[[i]][hd.list[[i]] <= tol] <- 0 # Set value less than tolerance to zero
  hd.list[[i]] <- as.dist(hd.list[[i]])
}

# Make a dataframe of horn distances
hd.frame <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  hd.frame[[i]] <- data.table(cbind(hd.list[[i]], env.dist))
}

hd.frame <- rbindlist(hd.frame, idcol = T)
setnames(hd.frame, c(".id", "V1"), c("sim", "horn.dist"))
hd.frame[, sim := as.factor(sim)]

# Calculate correlation 
hd.cor <- hd.frame[, .(spearman = cor(horn.dist, env.dist, method = "spearman"),
                       pearson = cor(horn.dist, env.dist)), by = sim]

# Summarize and write to file
# Spearman
hd.cor.spearman <- summary(hd.cor[, spearman])
hd.cor.spearman <- round(unclass(hd.cor.spearman), 5)
hd.cor.spearman <- data.table(t(hd.cor.spearman))

write.csv(hd.cor.spearman, 
          file = paste(here(), 
                       "/analysis/diffusion-map-high-turnover/",
                       "horn-spearman-summary.csv", 
                       sep = ""))

# Pearson
hd.cor.pearson <- summary(hd.cor[, pearson])
hd.cor.pearson <- unclass(hd.cor.pearson)
hd.cor.pearson <- data.table((t(hd.cor.pearson)))

write.csv(hd.cor.pearson, 
          file = paste(here(), 
                       "/analysis/diffusion-map-high-turnover/",
                       "horn-pearson-summary.csv", 
                       sep = ""))

# Plot diffusion distance vs environmental distance
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "horn-distance-plot.png", sep = ""),
    height = 2, width = 3, units = "in", res = 600)
ggplot(hd.frame, aes(x = env.dist, y = horn.dist)) + 
  geom_point(show.legend = F, size = 0.1) + 
  theme_bw() + 
  labs(x = "Environmental Distance", y = "Horn Distance")
dev.off()

#### Compare ordinations--------------------------------------------------------
# PCoA
pcoa.horn <- vector(mode = "list", length = n.sim)
for (i in 1:n.sim) {
  pcoa.horn[[i]] <- pcoa(hd.list[[i]])
}

# NMDS
nmds.horn <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  nmds.horn[[i]] <- metaMDS(hd.list[[i]], k = 2)
}

# Make reference matrix of environmental variables
env.mat <- as.matrix(cbind(grid.samp[1, ], grid.samp[2, ]))

# Procrustes statistic for horn distances
proc.pcoa <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  proc.pcoa[[i]] <- vegan::procrustes(env.mat, pcoa.horn[[i]]$vectors[, 1:2], 
                                      symmetric = T)
}

# Extract sum of squares and summarize
pcoa.ss <- vector(length = n.sim)

for (i in 1:n.sim) {
  pcoa.ss[i] <- proc.pcoa[[i]]$ss
}

# Summarize and write to file
pcoa.ss.summ <- summary(pcoa.ss)
pcoa.ss.summ <- round(unclass(pcoa.ss.summ), 5)
pcoa.ss.summ <- data.table(t(pcoa.ss.summ))

write.csv(pcoa.ss.summ, 
          file = paste(here(), 
                       "/analysis/diffusion-map-high-turnover/",
                       "pcoa-procrustes-summary.csv", 
                       sep = ""))

# Single plot with procrustes sum of squares error
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "pcoa-procrustes-plot.png", sep = ""),
    height = 3.5, width = 3.5, units = "in", res = 600)
plot(proc.pcoa[[1]], main = "")
dev.off()


# Procrustes statistic for NMDS
proc.nmds <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  proc.nmds[[i]] <- vegan::procrustes(env.mat, nmds.horn[[i]]$points, 
                                      symmetric = T)
}

# Extract sum of squares and summarize
nmds.ss <- vector(length = n.sim)

for (i in 1:n.sim) {
  nmds.ss[i] <- proc.nmds[[i]]$ss
}

# Summarize and write to file
nmds.ss.summ <- summary(nmds.ss)
nmds.ss.summ <- round(unclass(nmds.ss.summ), 5)
nmds.ss.summ <- data.table(t(nmds.ss.summ))

write.csv(nmds.ss.summ, 
          file = paste(here(), 
                       "/analysis/diffusion-map-high-turnover/",
                       "nmds-procrustes-summary.csv", 
                       sep = ""))

# Single plot with procrustes sum of squares error
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "nmds-procrustes-plot.png", sep = ""),
    height = 3.5, width = 3.5, units = "in", res = 600)
plot(proc.nmds[[1]], main = "")
dev.off()


# Procrustes statistic for diffusion maps
proc.diff <- vector(mode = "list", length = n.sim)

for (i in 1:n.sim) {
  proc.diff[[i]] <- vegan::procrustes(env.mat, 
                                      site.dim[sim == i, as.matrix(.(dim1, dim2))], 
                                      symmetric = T)
}

# Extract sum of squares and summarize
diff.ss <- vector(length = n.sim)

for (i in 1:n.sim) {
  diff.ss[i] <- proc.diff[[i]]$ss
}

# Summarize and write to file
diff.ss.summ <- summary(diff.ss)
diff.ss.summ <- round(unclass(diff.ss.summ), 5)
diff.ss.summ <- data.table(t(diff.ss.summ))

write.csv(diff.ss.summ, 
          file = paste(here(), 
                       "/analysis/diffusion-map-high-turnover/",
                       "diffmap-procrustes-summary.csv", 
                       sep = ""))

# Single plot with procrustes sum of squares error
png(filename = paste(here(), "/figures/diffusion-map-high-turnover/",
                     "diffmap-procrustes-plot.png", sep = ""),
    height = 3.5, width = 3.5, units = "in", res = 600)
plot(proc.diff[[1]], main = "")
dev.off()