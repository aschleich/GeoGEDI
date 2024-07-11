FROM rocker/geospatial:4.4
USER rstudio

WORKDIR /home/rstudio/scripts
COPY *.R .
CMD ["/bin/bash"]
