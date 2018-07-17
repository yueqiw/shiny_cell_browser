FROM r-base

RUN apt-get update && apt-get install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev libhdf5-dev

RUN mkdir -p /app/data

WORKDIR /app

VOLUME /app/data

RUN R -e 'install.packages(c("Seurat"))'
RUN R -e 'install.packages(c("rjson"))'
RUN R -e 'install.packages(c("shiny"))'
RUN R -e 'install.packages(c("tidyverse"))'

ADD . /app/

EXPOSE 4242

CMD R -e "shiny::runApp('./', host='0.0.0.0', port=4242)"
