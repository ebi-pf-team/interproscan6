# This image contains only dependencies for DeepTMHMM. It does not include a copy of DeepTMHMM
FROM python:3.9-slim

ENV TZ=Europe/London

RUN apt-get -y update \
 && apt-get -y install procps

RUN pip install --no-cache-dir cython==0.29.37 pkgconfig==1.5.5 \
 && pip install --no-cache-dir torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117

COPY requirements.txt /tmp/requirements.txt

RUN pip install --no-cache-dir -r /tmp/requirements.txt
