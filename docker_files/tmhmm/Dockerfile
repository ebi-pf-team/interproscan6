FROM python:3.7-slim
ENV TZ=Europe/London

RUN apt-get update && \
apt-get upgrade -y && \
apt-get install -y \
    procps
RUN pip install -U pybiolib
RUN biolib run DTU/DeepTMHMM --help