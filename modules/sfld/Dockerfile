FROM interpro/hmmer:3.3

LABEL name="interpro/sfld"
LABEL tag="ebiftp"

ENV PATH="/opt/sfld:${PATH}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        procps \
        tar && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/sfld
RUN curl -L -O https://ftp.ebi.ac.uk/pub/databases/interpro/databases/sfld/sfld_postprocess && \
    chmod +x ./sfld_postprocess
