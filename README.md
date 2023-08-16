# Dissimilarity Analysis Based on Diffusion Maps

This repository contains the code necessary to replicate the results from Gault et al. (2023) Dissimilarity analysis based on diffusion maps. 

### Running analysis in Docker container

We have created a Docker image that contains the same version of R (v4.1.2) and all packages used in the analysis. The image can be downloaded directly from Dockerhub or built from the Dockerfile. 

Pull the image directly from Dockerhub via the following:

`docker pull jagault/diffusion-map:rstudio4.1.2`

Or, from within the directory containing the Dockerfile, build the image via the following:

`docker build . -t jagault/diffusion-map:rstudio4.1.2`  

Once the image is built, it can be run within the main repo directory `/diffusion-map-manuscript` via the following command:

`docker run -d -p 8787:8787 -v /path/to/repo:/home/rstudio -e PASSWORD=yourpasswordhere jagault/diffusion-map:rstudio4.1.2`  

The use of `-d` makes the container run in the background. `-p 8787:8787` opens a port through which rstudio can be accessed through an internet browser. On a linux machine, rstudio can be accessed via `http://localhost:8787`.  

If you are using a Windows or Mac machine, you can access rstudio via `http://yourip:8787`. Your ip address should be displayed when opening up the Docker Quickstart terminal.  

Rstudio no longer supports the use of the default password when none is specified. The user must now specify a password following the `-e` tag. 

The use of `-v /path/to/repo:/home/rstudio` will set up a bind mount between the project repo and rstudio within the container. All files within the repo will be available within Rstudio and saved changes will persist outside of the container.