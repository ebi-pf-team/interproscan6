FROM ubuntu:latest
LABEL authors="Laise Florentino (lcf@ebi.ac.uk), Matthias Blum (mblum@ebi.ac.uk)"
ARG VERSION=6.0-95.0
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get update -y && \
    apt-get install -y wget python3.10 python3-pip unzip python3-requests

WORKDIR /opt/interproscan6
COPY subworkflows/ subworkflows/
COPY scripts/ scripts/
COPY modules/ modules/
COPY interproscan.nf interproscan.nf
COPY nextflow.config nextflow.config
COPY README.md README.md
