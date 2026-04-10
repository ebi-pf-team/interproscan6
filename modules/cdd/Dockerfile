FROM ubuntu:24.04

LABEL name="interpro/rpsproc"
LABEL tag="0.5.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/blast:/opt/rpsbproc/RpsbProc-x64-linux:${PATH}"
ENV TZ=Europe/London

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        libdw1 \
        libgomp1 \
        procps \
        tar && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/blast
RUN curl -L -O https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    tar -zxpf ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    mv ncbi-blast-2.15.0+/bin/rpsblast rpsblast && \
    rm -rf ncbi-blast-2.15.0+

WORKDIR /opt/rpsbproc
RUN curl -L -O https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz && \
    tar -xzf RpsbProc-x64-linux.tar.gz && \
    rm RpsbProc-x64-linux.tar.gz
