FROM --platform=linux/amd64 ubuntu:latest as interproscan-base
LABEL authors="Laise Florentino (lcf@ebi.ac.uk), Emma Hobbs (ehobbs@ebi.ac.uk), Matthias Blum (mblum@ebi.ac.uk)"
ARG VERSION=6.0-95.0
ENV TZ=Europe/London
ENV NXF_ANSI_LOG=false
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get update -y && \
    apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC && \
    apt-get install -y autoconf automake autotools-dev bison build-essential cmake curl \
    flex git libcurl3-gnutls libdivsufsort3 liblmdb0 libdw1 libgomp1 libnghttp2-dev \
    libssl-dev libtool nghttp2 procps python3.10 python3-venv python3-pip python3-requests \
    tar unzip zlib1g-dev

# Pull pftools for HAMAP and PROSITE
FROM sibswiss/pftools as pftools

# Final image with InterProScan and pftoools
FROM interproscan-base
COPY --from=pftools / /opt/pftools
ENV PATH="/opt/pftools/usr/local/bin:${PATH}"

# Install NCBI BLAST, only rpsblast (for CDD)
# Don't pull the NCBI BLAST image has its BIG - Just get the bits we need
WORKDIR /opt/blast
RUN curl -L -O https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    tar -zxpf ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.15.0+-x64-linux.tar.gz

# Install RpsbProc for CDD post-processing
WORKDIR /opt/rpsbproc
RUN curl -L -O https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz && \
    tar -xzf RpsbProc-x64-linux.tar.gz && \
    rm RpsbProc-x64-linux.tar.gz

# Install HMMER3
WORKDIR /opt/
RUN mkdir /opt/hmmer3 && \
    curl -L -O http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz && \
    tar -xzf hmmer-3.3.tar.gz && \
    rm hmmer-3.3.tar.gz && \
    cd hmmer-3.3 && \
    ./configure --prefix /opt/hmmer3 && \
    make && \
    make install

# Install HMMER2 (for SMART)
WORKDIR /opt/
RUN mkdir /opt/hmmer2 && \
    curl -L -O http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz && \
    tar -xzf hmmer-2.3.2.tar.gz && \
    rm hmmer-2.3.2.tar.gz && \
    cd hmmer-2.3.2 && \
    ./configure --prefix /opt/hmmer2 && \
    make && \
    make install && \
    rm -rf /opt/hmmer-2.3.2

# Install easel for predicting open reading frames (ORFs)
RUN cd /opt/hmmer-3.3/easel && \
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
RUN pip install biopython==1.83

# Install Cath-tools for Gene3D and FunFam
WORKDIR /opt/cath-tools
RUN curl -L -o cath-resolve-hits https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.10/cath-resolve-hits.ubuntu-20.04
RUN chmod +x cath-resolve-hits

WORKDIR /opt/interproscan6
COPY interproscan/ interproscan/
COPY utilities/ utilities/
COPY tests/ tests/
COPY interproscan.nf interproscan.nf
COPY nextflow.config nextflow.config
COPY README.md README.md
COPY LICENSE LICENSE
