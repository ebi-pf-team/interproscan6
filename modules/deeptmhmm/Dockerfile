# This image contains only dependencies for DeepTMHMM. It does not include a copy of DeepTMHMM
FROM python:3.9-slim

LABEL name="interpro/deeptmhmm"
LABEL tag="1.0"

ENV TZ=Europe/London

RUN apt-get -y update \
 && apt-get -y install procps \
 && apt-get clean

RUN pip install --no-cache-dir cython==0.29.37 pkgconfig==1.5.5 \
 && pip install --no-cache-dir torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117

RUN pip install --no-cache-dir \
    biopython==1.79 \
    certifi==2022.9.24 \
    charset-normalizer==2.1.1 \
    cycler==0.11.0 \
    fair-esm==0.4.0 \
    fonttools==4.38.0 \
    future==0.18.2 \
    h5py==3.7.0 \
    idna==3.15 \
    kiwisolver==1.4.4 \
    matplotlib==3.5.2 \
    numpy==1.23.4 \
    packaging==21.3 \
    PeptideBuilder==1.1.0 \
    Pillow==12.2.0 \
    pyparsing==3.0.9 \
    python-dateutil==2.8.2 \
    requests==2.33.0 \
    six==1.16.0 \
    tqdm==4.64.1 \
    urllib3==2.7.0
