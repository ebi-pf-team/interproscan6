FROM ubuntu:22.04
ENV TZ="UTC"
WORKDIR /cath-tools

RUN apt-get update && apt-get install -y wget
RUN wget -O cath-resolve-hits https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.10/cath-resolve-hits.ubuntu-20.04
RUN chmod +x cath-resolve-hits
