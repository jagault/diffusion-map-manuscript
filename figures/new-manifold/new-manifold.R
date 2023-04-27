#### Load required libraries and helper functions
library(here)
library(data.table)
library(vegan)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(gg3D)
library(viridis)

source(here("R", "diffusion_map.R"))
source(here("R", "horn_overlap.R"))
### Generate response curves----------------------------------------------------
set.seed(10)

M <- 3 # size of species inventory
N <- 50 # sample size
# nichepos <- seq(0.1, 0.9, length.out = M)
# nichewidth <- runif(M, 0.2, 0.3)

nichepos <- c(0.4, 0.5, 0.6)
nichewidth <- c(0.2, 0.4, 0.2)

# nichepos <- c(0.4, 0.5, 0.6)
# nichewidth <- c(0.15, 0.4, 0.15)



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
png(filename = paste(here(), "/figures/new-manifold/",
                     "response-curves.png", sep = ""),
    height = 4.5, width = 5.5, units = "in", res = 600)
matplot(as.data.frame(xseq), as.data.frame(niches), type='l', 
        xlab='Environmental Variable', ylab='Expected Species Density',
        col = color, lwd = 3)
dev.off()

#### Sample from response curves to generate assemblages------------------------
xsample <- seq(0, 1, length.out = N) # evenly spaced environmental samples
dens <- matrix(nrow = N, ncol = M)
abundance <-  matrix(nrow = N, ncol = M)
area <- 20

# Single realizations of relative fitness
# for (col in 1:M){
#   dens[,col] <- lambda(xsample,  nichepos[col], nichewidth[col])
#   abundance[,col] <- dens[,col]*area
#   # abundance[,col] <- rpois(N, dens[,col]*area)
# }
# 
# relabundance <- t(apply(abundance, 1, function(x) x/sum(x)))

# Multiple realizations for each environmental value
n.sim <- 1000 # Number of realizations

set.seed(550)
comm <- vector(mode = "list", length = n.sim)
for (i in 1:n.sim) {
  
  for (col in 1:M){
    dens[,col] <- lambda(xsample,  nichepos[col], nichewidth[col])
    # abundance[,col] <- round(dens[,col]*area)
    abundance[,col] <- rpois(N, dens[,col]*area)
  }
  comm[[i]] <- data.table(abundance)
}

comm <- rbindlist(comm, idcol = T)
comm[, env := rep(1:N, n.sim)]
comm <- melt(comm, id.vars = c(".id", "env"), variable.name = "spec", value.name = "abund")
comm <- comm[, sum(abund), by = c("spec", "env")]
comm <- dcast(comm, env ~ spec)
setnames(comm, c("V1", "V2", "V3"), c("A", "B", "C"))

# comm[, c("A", "B", "C") := data.table(t(apply(.SD, 1, function(x) x/sum(x)))), 
#      .SDcols = c("A", "B", "C")]
# relabundance <- comm

relabundance <- comm[, t(apply(.SD, 1, function(x) x/sum(x))), 
.SDcols = c("A", "B", "C")]


#### Plot assemblage compositions-----------------------------------------------
# Format data
cm <- relabundance
rownames(cm) <- paste("S", 1:nrow(cm), sep = "")
colnames(cm) <- LETTERS[1:ncol(cm)]

# Set as datatable for plotting
cm.dt <- data.table(cm, keep.rownames = T)
cm.dt[, env := xsample]
setnames(cm.dt, "rn", "site")

# Plot
png(filename = paste(here(), "/figures/new-manifold/",
                     "manifold.png", sep = ""),
    height = 7, width = 7, units = "in", res = 600)
ggplot(cm.dt, aes(x = A, y = B, z = C, color = env)) + 
  theme_void() + 
  scale_color_viridis(option = "plasma", name = "Env") + 
  axes_3D() + 
  stat_3D() + 
  labs_3D(labs = c("Species A", "Species B", " Species C"),
          hjust=c(0.0001, 1.2, 1), vjust=c(2, 2, -0.2), angle=c(45, -45, 0))
dev.off()


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
env <- xsample
site.dim[, env := env]
# site.dim[, richness := apply(cm > 0, sum, MARGIN = 1)]
setcolorder(site.dim, c("site", "dim1", "dim2", "dim3", "env"))

# Flip axes to for easier visualization
# site.dim[, dim1 := dim1*-1]
# site.dim[, dim2 := dim2*-1]

# Plot diffmap weighted by eigenvalues
png(filename = paste(here(), "/figures/new-manifold/",
                     "diffmap.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(site.dim, aes(x = dim1, y = dim2)) + 
  geom_point(aes(color = env), size = 2) + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  scale_color_viridis(option = "plasma", name = "Env. Var.") + 
  theme_bw() + 
  theme(legend.key.size = unit(0.2, "in"),
        text = element_text(size = 15)) +
  xlim(-110, 110) + ylim(-110, 110)
dev.off()

# Species Responses: Diffusion Map----------------------------------------------
# Format data
# a <- abundance
# rownames(a) <- paste("S", 1:nrow(a), sep = "")
# colnames(a) <- LETTERS[1:ncol(a)]

# Set as datatable for plotting
# a.dt <- data.table(a, keep.rownames = T)
# setnames(a.dt, "rn", "site")

a.dt <- comm[, .(A, B, C)]
a.dt[, site := paste("S", 1:nrow(a.dt), sep = "")]

setkey(a.dt, site)
setkey(site.dim, site)

res <- merge(site.dim, a.dt)
res.melt <- melt(res, measure.vars = c("A", "B", "C"), variable.name = "spec", 
                 value.name = "abund")

png(filename = paste(here(), "/figures/new-manifold/",
                     "estimated-response-curves.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(res.melt, aes(x = dim1, y = abund, color = spec)) + 
  geom_point() +
  geom_line() +
  scale_color_manual(name = "Species",
                     breaks = c("A", "B", "C"), 
                     values = brewer.pal(10, "Paired")[1:3]) +
  labs(x = "Dimension 1", y = "Abundance") +
  theme_bw() +
  theme(legend.key.size = unit(0.2, "in"),
        text = element_text(size = 15))
dev.off()

#### NMDS-----------------------------------------------------------------------
# Calculate Horn dissimilarity
hs <- HornMat(cm)
hd <- 1 - hs
tol <- 1e-10 # Set tolerance for converting small floats to 0
hd[hd <= tol] <- 0

nm <- metaMDS(hd, k = 2)

nm.dt <- data.table(nm$points)
nm.dt[, env := env]
nm.dt[, site := paste("S", 1:nrow(a.dt), sep = "")]

png(filename = paste(here(), "/figures/new-manifold/",
                     "nmds.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(nm.dt, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = env), size = 2) + 
  labs(x = "MDS1", y = "MDS2") + 
  scale_color_viridis(option = "plasma", name = "Env. Var.") + 
  theme_bw() + 
  theme(legend.key.size = unit(0.2, "in"),
        text = element_text(size = 15)) 
dev.off()

#### Species Response Curves: NMDS----------------------------------------------
# a.dt <- comm[, .(A, B, C)]
# a.dt[, site := paste("S", 1:nrow(a.dt), sep = "")]
# 
# setkey(a.dt, site)
# setkey(nm.dt, site)
# 
# res <- merge(nm.dt, a.dt)
# res.melt <- melt(res, measure.vars = c("A", "B", "C"), variable.name = "spec", 
#                  value.name = "abund")
# 
# ggplot(res.melt, aes(x = MDS1, y = abund, color = spec)) + 
#   geom_point() +
#   geom_line() +
#   scale_color_manual(breaks = c("A", "B", "C"), 
#                      values = brewer.pal(10, "Paired")[3:1]) +
#   theme_bw()

#### Distance-------------------------------------------------------------------
setkey(site.dim, env)
# Environmental distances
env.dist <- round(dist(site.dim[, .(env)]), 4)

# Diffusion distance
diff.dist <- round(dist(as.matrix(site.dim[, .(dim1)])), 4)
diff.env <- data.table(cbind(diff.dist, env.dist))

diff.cor <- diff.env[, .(spearman = cor(diff.dist, env.dist, 
                                        method = "spearman"),
                         pearson = cor(diff.dist, env.dist))]

png(filename = paste(here(), "/figures/new-manifold/",
                     "diffusion-distance-plot.png", sep = ""),
    height = 4, width = 6, units = "in", res = 600)
ggplot(diff.env, aes(x = env.dist, y = diff.dist)) + 
  geom_point(size = 0.5) + 
  theme_bw() +
  labs(x = "Environmental Distance", y = "Diffusion Distance") +
  theme(text = element_text(size = 15))
dev.off()

# # NMDS distance
# nmds.dist <- round(dist(nm$points), 4)
# nmds.env <- data.table(cbind(nmds.dist, env.dist))
# 
# nmds.cor <- nmds.env[, .(spearman = cor(nmds.dist, env.dist, 
#                                         method = "spearman"),
#                          pearson = cor(nmds.dist, env.dist))]
# 
# ggplot(nmds.env, aes(x = env.dist, y = nmds.dist)) + 
#   geom_point() + 
#   theme_bw() +
#   labs(x = "Environmental Distance", y = "NMDS Distance")
# 
# # Horn distance
# l.hd <- hd[lower.tri(hd, diag = F)]
# hd.env <- data.table(cbind(l.hd, env.dist))
# 
# ggplot(hd.env, aes(x = env.dist, y = nmds.dist)) + 
#   geom_point() + 
#   theme_bw() +
#   labs(x = "Environmental Distance", y = "Horn Distance")
