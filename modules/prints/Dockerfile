FROM ubuntu:24.04

LABEL name="interpro/fingerprintscan"
LABEL tag="3.597.ebiftp"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/fingerprintscan:${PATH}"
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

WORKDIR /opt/fingerprintscan
RUN curl -L -O https://ftp.ebi.ac.uk/pub/databases/interpro/databases/prints/fingerPRINTScan && \
    chmod +x ./fingerPRINTScan
