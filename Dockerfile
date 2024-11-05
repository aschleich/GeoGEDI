FROM rocker/geospatial:4.4

RUN apt-get update -y && apt-get install parallel -y --no-install-recommends

RUN R -q -e 'install.packages(c("whitebox", "ModelMetrics", "arrow", "foreach", "doParallel"))' \
    && strip /usr/local/lib/R/site-library/*/libs/*.so

USER ubuntu
WORKDIR /home/ubuntu
RUN R -q -e 'library(whitebox) ; install_whitebox()'

ENV GDAL_CACHEMAX=1024
ENV GDAL_DISABLE_READDIR_ON_OPEN=true

ENV PATH="/home/ubuntu/.local/bin:$PATH"
COPY *.R /usr/local/bin/
CMD ["/bin/bash"]
