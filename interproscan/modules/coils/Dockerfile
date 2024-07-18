FROM --platform=linux/amd64 ubuntu:latest as interproscan-base
LABEL authors="Laise Florentino (lcf@ebi.ac.uk), Emma Hobbs (ehobbs@ebi.ac.uk), Matthias Blum (mblum@ebi.ac.uk)"
ARG VERSION=6.0-95.0
ENV TZ=Europe/London
ENV NXF_ANSI_LOG=false
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get update -y && \
    apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC && \
    apt-get install -y autoconf automake autotools-dev bison build-essential cmake curl \
    flex git libcurl3-gnutls libdivsufsort3 liblmdb0 libdw1 libgomp1 libnghttp2-dev \
    libssl-dev libtool nghttp2 procps python3.10 python3-venv python3-pip python3-requests \
    tar unzip zlib1g-dev

WORKDIR /opt/coils
RUN curl -O -L http://www.russelllab.org/coils/coils.tar.gz && \
    tar -xzf coils.tar.gz &&\
    rm coils.tar.gz && \
    cd coils
