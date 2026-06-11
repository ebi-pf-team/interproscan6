FROM python:3.7-slim

LABEL name="interpro/signalp"
LABEL tag="6.0i"

ENV TZ=Europe/London
RUN apt-get update -y && \
    apt-get install -y procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    matplotlib~=3.5.3 \
    numpy~=1.21.6 \
    torch~=1.13.1 \
    tqdm~=4.67.1
