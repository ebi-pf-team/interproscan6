FROM --platform=linux/amd64 ubuntu:latest
ENV TZ=Europe/London
WORKDIR /opt/cdd

RUN apt-get update && apt-get install -y \
    curl \
    tar \
    libdw1 \
    libgomp1
RUN curl -O https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz && \
    tar -xzf RpsbProc-x64-linux.tar.gz && \
    rm RpsbProc-x64-linux.tar.gz
