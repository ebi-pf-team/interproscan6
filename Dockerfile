# Base image with dependencies
FROM --platform=linux/amd64 ubuntu:latest as interproscan-base
LABEL authors="Laise Florentino (lcf@ebi.ac.uk), Matthias Blum (mblum@ebi.ac.uk)"
ARG VERSION=6.0-95.0
ENV TZ=Europe/London
ENV NXF_ANSI_LOG=false
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get update -y && \
    apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC && \
    apt-get install -y autoconf automake autotools-dev bison build-essential cmake curl flex git libcurl3-gnutls libdivsufsort3 liblmdb0 libdw1 libgomp1 libnghttp2-dev libssl-dev libtool nghttp2 python3.10 python3-venv python3-pip python3-requests tar unzip wget zlib1g-dev

# Pull NCBI BLAST (for CDD)
FROM ncbi/blast as blast

# Pull HMMER container
FROM biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1 as hmmer

# Final image with InterProScan, BLAST, and HMMER
FROM interproscan-base
COPY --from=blast / /blast
COPY --from=hmmer / /hmmer

# Install easel for predicting open reading frames (ORFs)
WORKDIR /opt/easel
RUN git clone https://github.com/EddyRivasLab/easel && \
    cd easel && \
    autoconf && \
    ./configure && \
    make && \
    make check

# Install epa-ng and biopython for Panther post-processing
WORKDIR /opt/
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN git clone https://github.com/pierrebarbera/epa-ng
RUN cd epa-ng && make
RUN pip install biopython

# Install RpsbProc for CDD post-processing
RUN curl -O https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz && \
    tar -xzf RpsbProc-x64-linux.tar.gz && \
    rm RpsbProc-x64-linux.tar.gz

WORKDIR /opt/interproscan6
COPY . .
