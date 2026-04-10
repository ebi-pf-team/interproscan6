FROM interpro/hmmer:3.3

LABEL name="interpro/cath"
LABEL tag="0.16.10"

ENV PATH="/opt/cath-tools:${PATH}"

WORKDIR /opt/cath-tools

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN curl -L -o /opt/cath-tools/cath-resolve-hits https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.10/cath-resolve-hits.ubuntu-20.04 && \
    chmod +x /opt/cath-tools/cath-resolve-hits
