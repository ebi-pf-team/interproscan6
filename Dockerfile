FROM ubuntu:22.04
LABEL authors="Laise Florentino (lcf@ebi.ac.uk), Matthias Blum (mblum@ebi.ac.uk)"
ARG VERSION=6.0-95.0
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get update -y && \
    apt-get install -y wget python3.10 python3-pip hmmer unzip

RUN python3 -m pip install --upgrade pip && \
    pip install biopython pyhmmer

WORKDIR /opt/interproscan6
COPY . .
ENTRYPOINT ["python3", "scripts/members/cdd.py"]
