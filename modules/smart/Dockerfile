FROM interpro/hmmer:3.3

LABEL name="interpro/smart"
LABEL tag="1.0"

ENV PATH="/opt/hmmer2/bin:${PATH}"

WORKDIR /opt/hmmer2
RUN curl -L -O http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz && \
    tar -xzf hmmer-2.3.2.tar.gz && \
    rm hmmer-2.3.2.tar.gz && \
    cd hmmer-2.3.2 && \
    ./configure --enable-threads --prefix /opt/hmmer2 && \
    make && \
    make install && \
    cd /opt/hmmer2 && \
    rm -rf hmmer-2.3.2 hmmer2/man man  && \
    find bin -type f ! -name 'hmmpfam' -delete
