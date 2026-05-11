FROM ubuntu:24.04

LABEL name="interpro/ncoils"
LABEL tag="2.2.1"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/ncoils:${PATH}"
ENV TZ=Europe/London

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        procps \
        tar && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/ncoils
RUN curl -L -O https://github.com/ebi-pf-team/ncoils/archive/refs/tags/2.2.1.tar.gz && \
    tar -zxpf 2.2.1.tar.gz && \
    cd ncoils-2.2.1 && \
    make && \
    chmod +x ./ncoils && \
    mv ./ncoils .. && \
    cd .. && \
    rm -rf ncoils-2.2.1
