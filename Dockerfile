FROM rocker/verse@sha256:21bb3cf9b3a843f590c4a4e2f916d4da42e3c231db31cfb23516e0ac54dba91a

ARG CONTAINER_USER=rstudio
ARG R_PACKAGE=CFDscoRe
RUN apt update && apt-get install -y g++-9
COPY . /build/${R_PACKAGE}
RUN cd /build/${R_PACKAGE} && R --quiet -e "devtools::document(); devtools::install()"
USER ${CONTAINER_USER}
