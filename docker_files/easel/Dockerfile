FROM ubuntu:latest
ENV TZ=Europe/London
WORKDIR /opt/easel
RUN apt-get update
RUN apt-get install -y curl build-essential nghttp2 libnghttp2-dev libssl-dev libcurl3-gnutls
RUN apt-get install -y git
RUN apt-get install -y autoconf


RUN git clone https://github.com/EddyRivasLab/easel && \
cd easel && \
autoconf && \
./configure && \
make && \
make check