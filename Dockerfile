FROM rocker/geospatial:4.4

RUN apt-get update -y && apt-get install parallel -y --no-install-recommends

RUN useradd -u 1001 -s /bin/bash -m geogedi
USER geogedi
WORKDIR /home/geogedi/scripts

COPY *.R .
CMD ["/bin/bash"]
