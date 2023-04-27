#### Load required libraries and helper functions
library(here)
library(ggplot2)
library(data.table)
library(vegan)


source(here("R", "diffusion_map.R"))
source(here("R", "horn_overlap.R"))

#### Load data------------------------------------------------------------------
pyke <- fread(file = paste(here(), "/analysis/empirical-example/pyke-2001.csv", 
                           sep = ""), 
              header = T)

con.site <- fread(file = paste(here(), 
                               "/analysis/empirical-example/condit-site-info.csv", 
                               sep = ""), 
                  header = T, check.names = T)

con.spec <- fread(file = paste(here(), 
                               "/analysis/empirical-example/condit-species-matrix.csv", 
                               sep = ""), 
                  header = T)

### Merge site info-------------------------------------------------------------
pyke[, coord := paste(pyke$utm.x, pyke$utm.y, sep = "-")]

con.site[, V8 := NULL]
con.site[, EW.coord := round(con.site$EW.coord)]
con.site[, NS.coord := round(con.site$NS.coord)]
con.site[, coord := paste(con.site$EW.coord, con.site$NS.coord, sep = "-")]

setkey(pyke, coord)
setkey(con.site, coord)

sites <- merge(pyke, con.site)
sites

### Select BCI sites to drop----------------------------------------------------
bs <- paste(rep("B", 50), c(0:49), sep = "")
con.site[site.no. %in% bs, ]

ggplot(data = con.site[site.no. %in% bs, ], aes(x = EW.coord, y = NS.coord)) +
  geom_text(aes(label = site.no.)) + 
  theme_bw()
# They kept the corners so keep sites 0, 4, 45, and 49


# Diffusion map-----------------------------------------------------------------
con.spec
cdat <- copy(con.spec)
drop <- paste(rep("BCI", 50), c(0:49), sep = "")
drop <- drop[!(drop %in% c("BCI0", "BCI4", "BCI45", "BCI49"))]
cdat[, c(drop) := NULL]

cdat[, species := NULL]
cm <- as.matrix(cdat)
rownames(cm) <- con.spec[, species]
cm <- t(cm)

cm <- decostand(cm, method = "total")

### Diffusion map
# Calculate dissimilarity between samples
d <- vegdist(cm, method = "horn")
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

# Add site column
site.dim[, site := rownames(cm)]

# Plot
ggplot(site.dim, aes(x = dim1, y = dim2)) +
  geom_point() + 
  labs(x = "Dimension 1", y = "Dimension 2") + 
  theme_bw() +
  xlim(-2.5, 3.5) + ylim(-2.5, 3.5)

ggplot(site.dim, aes(x = dim1, y = dim2)) +
  geom_text(aes(label = site)) + 
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

# Plot geo coords
ggplot(data = con.site, aes(x = EW.coord, y = NS.coord)) +
  geom_text(aes(label = site.no.)) + 
  theme_bw()

bs <- paste(rep("B", 50), c(0:49), sep = "")
con.site[site.no. %in% bs, ]

ggplot(data = con.site[site.no. %in% bs, ], aes(x = EW.coord, y = NS.coord)) +
  geom_text(aes(label = site.no.)) + 
  theme_bw()
