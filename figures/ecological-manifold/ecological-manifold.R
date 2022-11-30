#### Load required libraries and helper functions
library(here)
library(data.table)
library(vegan)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(gg3D)

### Generate response curves----------------------------------------------------
set.seed(10)

M <- 3 # size of species inventory
N <- 25 # sample size
nichepos <- seq(0.1, 0.9, length.out = M)
nichewidth <- runif(M, 0.2, 0.3)

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
png(filename = paste(here(), "/figures/ecological-manifold/",
                     "repsonse-curves.png", sep = ""),
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
for (col in 1:M){
  dens[,col] <- lambda(xsample,  nichepos[col], nichewidth[col])
  # abundance[,col] <- round(dens[,col]*area)
  abundance[,col] <- rpois(N, dens[,col]*area)
}

relabundance <- t(apply(abundance, 1, function(x) x/sum(x)))

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
png(filename = paste(here(), "/figures/ecological-manifold/",
                     "manifold.png", sep = ""),
    height = 7, width = 7, units = "in", res = 600)
ggplot(cm.dt, aes(x = A, y = B, z = C, color = env)) + 
  theme_void() + 
  scale_color_continuous(type = "viridis", name = "Env") + 
  axes_3D() + 
  stat_3D() + 
  labs_3D(labs = c("Species A", "Species B", " Species C"),
          hjust=c(0.4,0.8,1), vjust=c(2, 2.1, -0.2), angle=c(0, 0, 90))
dev.off()