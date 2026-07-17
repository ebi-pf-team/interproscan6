FROM ubuntu:26.04 AS builder

LABEL name="interpro/nail"
LABEL tag="0.6.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# nail 0.6.0 is tested against MMseqs2 Release 15-6f452
ARG NAIL_VERSION=nail-v0.6.0
ARG MMSEQS_VERSION=15-6f452
# SIMD instruction set used to build MMseqs2 (AVX2 for modern x86-64;
# override with SSE4_1 for older CPUs).
ARG MMSEQS_SIMD=AVX2

ENV PATH="/root/.cargo/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        libbz2-dev \
        xxd \
        zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# MMseqs2 (nail's alignment prefilter)
# Built from the release tarball (VERSION_OVERRIDE replaces the git-derived
# version string, so no git checkout is needed here).
RUN mkdir -p /opt/src && \
    curl -L -o /tmp/mmseqs.tar.gz \
        https://github.com/soedinglab/MMseqs2/archive/refs/tags/${MMSEQS_VERSION}.tar.gz && \
    tar -xzf /tmp/mmseqs.tar.gz -C /opt/src && \
    rm /tmp/mmseqs.tar.gz && \
    cmake -S /opt/src/MMseqs2-${MMSEQS_VERSION} -B /opt/src/MMseqs2-${MMSEQS_VERSION}/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DHAVE_${MMSEQS_SIMD}=1 \
        -DVERSION_OVERRIDE=${MMSEQS_VERSION} \
        -DCMAKE_INSTALL_PREFIX=/opt/mmseqs && \
    cmake --build /opt/src/MMseqs2-${MMSEQS_VERSION}/build --target install -j"$(nproc)" && \
    rm -rf /opt/src/MMseqs2-${MMSEQS_VERSION}

# Rust toolchain + nail
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \
        sh -s -- -y --profile minimal --default-toolchain stable

RUN mkdir -p /opt/src && \
    curl -L -o /tmp/nail.tar.gz \
        https://github.com/TravisWheelerLab/nail/archive/refs/tags/${NAIL_VERSION}.tar.gz && \
    tar -xzf /tmp/nail.tar.gz -C /opt/src && \
    rm /tmp/nail.tar.gz && \
    cargo build --release --manifest-path /opt/src/nail-${NAIL_VERSION}/Cargo.toml && \
    mkdir -p /opt/nail/bin && \
    cp /opt/src/nail-${NAIL_VERSION}/target/release/nail /opt/nail/bin/nail && \
    rm -rf /opt/src/nail-${NAIL_VERSION} /root/.cargo/registry

# Runtime stage: keep only the compiled binaries and their runtime deps
FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London
ENV PATH="/opt/nail/bin:/opt/mmseqs/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libatomic1 \
        libbz2-1.0 \
        libgomp1 \
        procps \
        zlib1g && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/nail/bin/nail /opt/nail/bin/nail
COPY --from=builder /opt/mmseqs/bin/mmseqs /opt/mmseqs/bin/mmseqs
