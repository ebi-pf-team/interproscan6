FROM debian:bookworm-slim

LABEL name="interpro/phobius"
LABEL tag="1.01"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        perl \
        procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work