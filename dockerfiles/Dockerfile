FROM rocker/rstudio:4.1.2

RUN apt-get update && apt-get install -y \
    r-cran-tcltk2

RUN install2.r \
    data.table \
    here \
    vegan \
    ape \
    ggplot2 \
    RColorBrewer \
    devtools \
    pracma

RUN R -e "devtools::install_github('AckerDWM/gg3D')"
