FROM ubuntu:24.04

LABEL name="interpro/hmmer"
LABEL tag="3.3"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/hmmer3/bin:${PATH}"
ENV TZ=Europe/London

WORKDIR /opt

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        procps \
        tar && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/hmmer3 && \
    curl -L -O http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz && \
    tar -xzf hmmer-3.3.tar.gz && \
    rm hmmer-3.3.tar.gz && \
    cd hmmer-3.3 && \
    ./configure --prefix /opt/hmmer3 && \
    make -j"$(nproc)" && \
    make install && \
    cd /opt && \
    rm -rf hmmer-3.3 /opt/hmmer3/share
