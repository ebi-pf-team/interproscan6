FROM interpro/hmmer:3.3

LABEL name="interpro/epa-ng"
LABEL tag="0.3.8"

ENV PATH="/opt/epa-ng:${PATH}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        autotools-dev \
        bison \
        cmake \
        flex \
        libtool \
        zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/epa-ng
RUN curl -L -O https://github.com/pierrebarbera/epa-ng/archive/refs/tags/v0.3.8.tar.gz && \
    tar -xf v0.3.8.tar.gz && \
    cd epa-ng-0.3.8 && \
    make -j"$(nproc)" && \
    mv bin/epa-ng /opt/epa-ng/ && \
    cd /opt/epa-ng && \
    rm -rf epa-ng-0.3.8/ v0.3.8.tar.gz
