FROM --platform=linux/amd64 ubuntu:latest as interproscan-base
LABEL authors="Laise Florentino (lcf@ebi.ac.uk), Emma Hobbs (ehobbs@ebi.ac.uk), Matthias Blum (mblum@ebi.ac.uk)"
ARG VERSION=6.0-95.0
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get update -y && \
    apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC && \
    apt-get install --no-install-recommends -y \
    autoconf automake autotools-dev bison build-essential cmake curl \
    flex git libcurl3-gnutls libdivsufsort3 liblmdb0 libdw1 libgomp1 libnghttp2-dev \
    libssl-dev libtool nghttp2 procps python3.10 python3-venv python3-pip \
    tar zlib1g-dev
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Pull pftools for HAMAP and PROSITE
FROM sibswiss/pftools as pftools

# Final image with InterProScan and pftools
FROM interproscan-base
COPY --from=pftools /var/lib/pftools/bin/pfsearchV3 /opt/pftools/pfsearchV3
COPY --from=pftools /var/lib/pftools/bin/pfscanV3 /opt/pftools/pfscanV3
COPY --from=pftools /var/lib/pftools/bin/ps_scan.pl /opt/pftools/ps_scan.pl

# Get binary for coils, only used by Coils
WORKDIR /opt/coils
RUN coils_url="https://raw.githubusercontent.com/ebi-pf-team/interproscan/master/core/jms-implementation/support-mini-x86-32/src/coils/ncoils/2.2.1" && \
    coils_files=$(curl -s https://api.github.com/repos/ebi-pf-team/interproscan/contents/core/jms-implementation/support-mini-x86-32/src/coils/ncoils/2.2.1 | grep "name" | cut -d '"' -f 4) && \
    for cfile in $coils_files; \
    do \
        curl -L -O "$coils_url/$cfile"; \
    done
RUN make
RUN find . -type f ! -name 'ncoils' -delete

# Install NCBI BLAST, only rpsblast (for CDD)
# Don't pull the NCBI BLAST image has its BIG - Just get the bits we need
WORKDIR /opt/blast
RUN curl -L -O https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    tar -zxpf ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN rm ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN mv ncbi-blast-2.15.0+/bin/rpsblast rpsblast && \
    rm -rf ncbi-blast-2.15.0+

# Install RpsbProc for CDD post-processing
WORKDIR /opt/rpsbproc
RUN curl -L -O https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz && \
    tar -xzf RpsbProc-x64-linux.tar.gz
RUN rm RpsbProc-x64-linux.tar.gz

# Install HMMER3
WORKDIR /opt/
RUN mkdir /opt/hmmer3 && \
    curl -L -O http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz && \
    tar -xzf hmmer-3.3.tar.gz
RUN rm hmmer-3.3.tar.gz && \
    cd hmmer-3.3 && \
    ./configure --prefix /opt/hmmer3 && \
    make && \
    make install
RUN rm -rf /opt/hmmer3/share
RUN cd /opt/hmmer3/bin && \
    find . -type f ! -name 'hmmpress' ! -name 'hmmsearch' ! -name 'hmmscan' -delete

# Install easel for predicting open reading frames (ORFs)
RUN mv /opt/hmmer-3.3/easel /opt/easel
WORKDIR /opt/easel
RUN autoconf && \
    ./configure && \
    make && \
    make check
RUN find . -mindepth 1 -maxdepth 1 ! -name 'miniapps' -exec rm -rf {} + && \
    cd miniapps && \
    find . -type f ! -name 'esl-translate*' -delete
RUN rm -rf /opt/hmmer-3.3

# Install HMMER2 (for SMART)
WORKDIR /opt/
RUN mkdir /opt/hmmer2 && \
    curl -L -O http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz && \
    tar -xzf hmmer-2.3.2.tar.gz
RUN rm hmmer-2.3.2.tar.gz && \
    cd hmmer-2.3.2 && \
    ./configure --prefix /opt/hmmer2 && \
    make && \
    make install
RUN rm -rf /opt/hmmer-2.3.2 && rm -rf /opt/hmmer2/man
RUN cd /opt/hmmer2/bin && \
    find . -type f ! -name 'hmmpfam' -delete

# Install epa-ng and biopython for Panther post-processing
WORKDIR /opt/
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN git clone https://github.com/pierrebarbera/epa-ng
RUN cd epa-ng && make && \
    rm Dockerfile
RUN pip install biopython==1.83

# Install Cath-tools for Gene3D and FunFam
WORKDIR /opt/cath-tools
RUN curl -L -o cath-resolve-hits https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.10/cath-resolve-hits.ubuntu-20.04
RUN chmod +x cath-resolve-hits

WORKDIR /opt/interproscan6
COPY interproscan/ interproscan/
COPY utilities/linter utilities/linter
COPY utilities/install_nf-test.sh utilities/install_nf-test.sh
COPY utilities/requirements-dev.txt utilities/requirements-dev.txt
COPY tests/ tests/
COPY interproscan.nf interproscan.nf
COPY nextflow.config nextflow.config
COPY README.md README.md
COPY LICENSE LICENSE

# remove tools not needed for running IPS6 to reduce the final image size
RUN apt-get remove -y autoconf automake autotools-dev bison build-essential flex git libcurl3-gnutls libdivsufsort3 liblmdb0 libdw1 libnghttp2-dev libssl-dev libtool nghttp2 zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
