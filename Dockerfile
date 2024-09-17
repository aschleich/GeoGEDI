FROM rocker/geospatial:4.4

RUN apt-get update -y && apt-get install parallel -y --no-install-recommends

RUN R -q -e 'install.packages(c("whitebox", "ModelMetrics", "arrow"))' \
    && strip /usr/local/lib/R/site-library/*/libs/*.so

RUN useradd -u 1001 -s /bin/bash -m geogedi
USER geogedi
RUN R -q -e 'library(whitebox) ; install_whitebox()'

COPY *.R /usr/local/bin/
CMD ["/bin/bash"]
