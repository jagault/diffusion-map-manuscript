# Save first laplacian to file
test.lap.remote <- site.lap[[1]]$laplacian
saveRDS(test.lap.remote, file = paste(here(), "/docker-test",
                              "/test-lap-remote.rds", sep = ""))

# Read in laplacian from local machine
test.lap.local <- readRDS(paste(here(), "/docker-test",
                                "/test-lap-local.rds", sep = ""))

# Eigendecomposition of both
test.eig.remote <- eigen(test.lap.remote)
test.eig.local <- eigen(test.lap.local)

# Test if they are equal
test.eig.remote$vectors
test.eig.local$vectors
all.equal(test.eig.remote$vectors, test.eig.local$vectors)

# Save site.dim for first sim
site.dim1.remote <- site.dim[sim == 1, ]
saveRDS(site.dim1.remote, file = paste(here(), "/docker-test",
                                       "/site-dim1-remote.rds", sep = ""))

# Read in local site dim
site.dim1.local <- readRDS(paste(here(), "/docker-test",
                                 "/site-dim1-local.rds", sep = ""))

# Test if they are equal
all.equal(site.dim1.local, site.dim1.remote)
head(site.dim1.local)
head(site.dim1.remote)

site.dim1.local[, dim1 := dim1*-1]
site.dim1.local[, dim2 := dim2*-1]

all.equal(site.dim1.local, site.dim1.remote)
head(site.dim1.local)
head(site.dim1.remote)
