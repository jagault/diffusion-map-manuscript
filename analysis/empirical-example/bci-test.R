#### Load required libraries and helper functions
library(here)
library(ggplot2)
library(data.table)
library(vegan)


source(here("R", "diffusion_map.R"))
source(here("R", "horn_overlap.R"))

#### Diffusion map--------------------------------------------------------------
# Load data
data(BCI.env)
data(BCI)
bci <- BCI

# Format data
bci <- decostand(bci, method = "total")

# Calculate dissimilarity between samples
d <- vegdist(bci, method = "horn")
d <- as.matrix(d)
s <- 1 - d # Convert to similarity for diffusion map function
diag(s) <- 0

# Get site laplacian and eigenvectors/values
site.lap <- Lap(sim = s, cmax = 10)
site.eig <- eigen(site.lap$laplacian)

# Weight eigenvectors by their respective eigenvalues
for (i in 1:ncol(site.eig$vectors)){
  site.eig$vectors[, i] <- site.eig$vectors[, i]/site.eig$values[i]
}

# Select number of dimensions to keep
k <- 1:10 
# Create vector of column names
cnames <- paste("dim", rev(k), sep = "")

# Select eigenvectors corresponding to the first k dimensions
site.dim <- data.table(site.eig$vectors[, rank(site.eig$values, 
                                               ties.method = "first") %in% 
                                          (k + 1)])

# Set column names and column order
setnames(site.dim, cnames)
setcolorder(site.dim, rev(cnames))

# Plot
ggplot(site.dim, aes(x = dim1, y = dim2)) +
  geom_point() + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  theme_bw()

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

