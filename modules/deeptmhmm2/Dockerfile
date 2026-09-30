FROM python:3.12-slim

LABEL name="interpro/deeptmhmm2"
LABEL tag="0.1.0-ffccb76"

WORKDIR /opt/deeptmhmm2

ENV TZ=Europe/London
ENV TORCH_HOME=/opt/torch
ENV TORCHINDUCTOR_CACHE_DIR=/tmp/torchinductor

RUN apt-get -y update \
 && apt-get -y install git procps \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* \
 && pip install --no-cache-dir torch==2.10.0 --index-url https://download.pytorch.org/whl/cu128 \
 && git clone --single-branch --branch add-threads-option https://github.com/matthiasblum/DeepTMHMM2.git \
 && cd DeepTMHMM2 \
 && git checkout ffccb760dbfd25af978cf3ee0b4955b9f590d9b6 \
 && pip install --no-cache-dir . \
 && cd .. \
 && mv DeepTMHMM2/checkpoints/* . \
 && dtm2 --model-dir /opt/deeptmhmm2 --device cpu --threads 2 --simplify-io DeepTMHMM2/sample.fasta /tmp/dtm2-test \
 && rm -rf DeepTMHMM2 /tmp/dtm2-test \
 && chmod -R a+rX /opt/torch /opt/deeptmhmm2
