FROM ubuntu:24.04

LABEL name="interpro/esl-translate"
LABEL tag="0.46"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/easel:${PATH}"
ENV TZ=Europe/London

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf \
        build-essential \
        ca-certificates \
        curl \
        procps \
        tar && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/easel
RUN curl -L -O http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz && \
    tar -xzf hmmer-3.3.tar.gz && \
    rm hmmer-3.3.tar.gz && \
    cd hmmer-3.3/easel && \
    autoconf && \
    ./configure && \
    make -j"$(nproc)" && \
    mv miniapps/esl-translate /opt/easel && \
    cd /opt/easel && \
    rm -rf hmmer-3.3 
