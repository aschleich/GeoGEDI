FROM rocker/geospatial:4.4

RUN apt-get update -y && apt-get install parallel -y --no-install-recommends

USER rstudio
WORKDIR /home/rstudio/scripts

COPY *.R .
CMD ["/bin/bash"]
