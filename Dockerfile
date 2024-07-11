FROM rocker/geospatial:4.4

RUN apt-get update -y && apt-get install parallel -y --no-install-recommends

RUN useradd -m -s -u 1001 -g 1001 /bin/bash geogedi
USER geogedi
WORKDIR /home/geogedi/scripts

COPY *.R .
CMD ["/bin/bash"]
