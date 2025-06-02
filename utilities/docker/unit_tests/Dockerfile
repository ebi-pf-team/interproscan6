FROM --platform=linux/amd64 ubuntu:latest
LABEL authors="Emma Hobbs (ehobbs@ebi.ac.uk)"
ARG VERSION=6.0-95.0
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ /etc/timezone && \
    apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y python3.10 python3-venv python3-pip

# make dirs to store the scripts to mock the bins used in the IRL image
WORKDIR /opt/

# Copy the files that will form the outputs/inputs to be used in the unit tests
WORKDIR /opt/interproscan6
COPY lib/ lib/
COPY modules/ modules/
COPY subworkflows/ subworkflows/
COPY utilities/ utilities/
COPY tests/ tests/
COPY main.nf main.nf
COPY nextflow.config nextflow.config
COPY README.md README.md
COPY LICENSE LICENSE
