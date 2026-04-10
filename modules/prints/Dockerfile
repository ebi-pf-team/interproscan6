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

# RUN curl -L -O <url> && \
#     tar -xf 3.597.tar.gz && \
#     cd fingerprintscan-3.597/src && \
#     chmod +x ./configure && \
#     ./configure && \
#     make && \
#     chmod +x ./fingerPRINTScan && \
#     mv ./fingerPRINTScan /opt/fingerprintscan && \
#     cd ../ && \
#     rm -rf ncoils-2.2.1
