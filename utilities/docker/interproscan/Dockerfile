FROM ubuntu:latest AS interproscan-base
LABEL authors="Laise Florentino (lcf@ebi.ac.uk), Emma Hobbs (ehobbs@ebi.ac.uk), Matthias Blum (mblum@ebi.ac.uk)"
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get update -y && \
    apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC && \
    apt-get install --no-install-recommends -y \
    autoconf automake autotools-dev bison build-essential cmake curl \
    flex git libcurl3-gnutls libdivsufsort3 liblmdb0 libdw1 libgomp1 libnghttp2-dev \
    libssl-dev libtool nghttp2 procps python3.10 python3-venv python3-pip python3-requests \
    tar zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# Pull pftools for HAMAP and PROSITE
FROM sibswiss/pftools AS pftools

# Final image with InterProScan and pftools
FROM interproscan-base
COPY --from=pftools /var/lib/pftools/bin/pfsearchV3 /opt/pftools/pfsearchV3
COPY --from=pftools /var/lib/pftools/bin/pfscanV3 /opt/pftools/pfscanV3

# Get binary for coils, only used by Coils
WORKDIR /opt/coils
RUN coils_url="https://raw.githubusercontent.com/ebi-pf-team/interproscan/master/core/jms-implementation/support-mini-x86-32/src/coils/ncoils/2.2.1" && \
    coils_files=$(curl -s https://api.github.com/repos/ebi-pf-team/interproscan/contents/core/jms-implementation/support-mini-x86-32/src/coils/ncoils/2.2.1 | grep "name" | cut -d '"' -f 4) && \
    for cfile in $coils_files; \
    do \
        curl -L -O "$coils_url/$cfile"; \
    done && \
    make && \
    find . -type f ! -name 'ncoils' -delete

# Install NCBI BLAST, only rpsblast (for CDD)
# Don't pull the NCBI BLAST image as its BIG! - Just get the bits we need
WORKDIR /opt/blast
RUN curl -L -O https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    tar -zxpf ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    mv ncbi-blast-2.15.0+/bin/rpsblast rpsblast && \
    rm -rf ncbi-blast-2.15.0+

# Install RpsbProc for CDD post-processing
WORKDIR /opt/rpsbproc
RUN curl -L -O https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz && \
    tar -xzf RpsbProc-x64-linux.tar.gz && \
    rm RpsbProc-x64-linux.tar.gz

# Install HMMER3 and easel (which is packaged into hmmer3)
# easel is used for predicting open reading frames (ORFs)
WORKDIR /opt/
RUN mkdir /opt/hmmer3 && \
    curl -L -O http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz && \
    tar -xzf hmmer-3.3.tar.gz && \
    rm hmmer-3.3.tar.gz && \
    cd hmmer-3.3 && \
    ./configure --prefix /opt/hmmer3 && \
    make && \
    make install && \
    rm -rf /opt/hmmer3/share && \
    cd /opt/hmmer3/bin && \
    find . -type f ! -name 'hmmpress' ! -name 'hmmsearch' ! -name 'hmmscan' -delete && \
    mv /opt/hmmer-3.3/easel /opt/easel && \
    cd /opt/easel && \
    autoconf && \
    ./configure && \
    make && \
    make check && \
    find . -mindepth 1 -maxdepth 1 ! -name 'miniapps' -exec rm -rf {} + && \
    cd miniapps && \
    find . -type f ! -name 'esl-translate*' -delete && \
    rm -rf /opt/hmmer-3.3

# Install HMMER2 (for SMART)
RUN mkdir /opt/hmmer2 && \
    curl -L -O http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz && \
    tar -xzf hmmer-2.3.2.tar.gz && \
    rm hmmer-2.3.2.tar.gz && \
    cd hmmer-2.3.2 && \
    ./configure --enable-threads --prefix /opt/hmmer2 && \
    make && \
    make install && \
    rm -rf /opt/hmmer-2.3.2 && \
    rm -rf /opt/hmmer2/man && \
    cd /opt/hmmer2/bin && \
    find . -type f ! -name 'hmmpfam' -delete

# Install epa-ng and biopython for Panther post-processing
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN git clone https://github.com/pierrebarbera/epa-ng && \
    cd epa-ng && make && \
    rm Dockerfile
RUN pip install biopython==1.83

# Install Cath-tools for Gene3D and FunFam
WORKDIR /opt/cath-tools
RUN curl -L -o cath-resolve-hits https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.10/cath-resolve-hits.ubuntu-20.04 && \
    chmod +x cath-resolve-hits

WORKDIR /opt/interproscan6

# Add all tool directories to PATH
ENV PATH="/opt/hmmer3/bin:/opt/cath-tools:/opt/blast:/opt/rpsbproc/RpsbProc-x64-linux:/opt/coils:/opt/easel/miniapps:/opt/pftools:/opt/epa-ng/bin:/opt/hmmer2/bin:$PATH"

# remove tools not needed for running IPS6 to reduce the final image size
RUN apt-get remove -y autoconf automake autotools-dev bison build-essential flex git libcurl3-gnutls libdivsufsort3 liblmdb0 libdw1 libnghttp2-dev libssl-dev libtool nghttp2 zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
