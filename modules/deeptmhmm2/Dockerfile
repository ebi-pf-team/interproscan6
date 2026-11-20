FROM python:3.12-slim

LABEL name="interpro/deeptmhmm2"
LABEL tag="0.1.0-526f6c6"

WORKDIR /opt/deeptmhmm2

ENV TZ=Europe/London
ENV TORCH_HOME=/opt/torch
ENV TORCHINDUCTOR_CACHE_DIR=/tmp/torchinductor

RUN apt-get -y update \
 && apt-get -y install git procps \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* \
 && pip install --no-cache-dir torch==2.10.0 --index-url https://download.pytorch.org/whl/cu128 \
 && git clone https://github.com/fteufel/DeepTMHMM2.git \
 && cd DeepTMHMM2 \
 && git checkout 526f6c60f93c0868c7315413da8e669d12373cc3 \
 && pip install --no-cache-dir . \
 && cd .. \
 && mv DeepTMHMM2/checkpoints/* . \
 && dtm2 --model-dir /opt/deeptmhmm2 --device cpu --threads 2 --simplify-io DeepTMHMM2/sample.fasta /tmp/dtm2-test \
 && rm -rf DeepTMHMM2 /tmp/dtm2-test \
 && chmod -R a+rX /opt/torch /opt/deeptmhmm2
