FROM --platform=linux/amd64 ubuntu:latest
LABEL authors="Emma Hobbs (ehobbs@ebi.ac.uk)"
ARG VERSION=6.0-95.0
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ /etc/timezone && \
    apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y python3.10 python3-venv python3-pip

# make dirs to store the scripts to mock the bins used in the IRL image
WORKDIR /opt/
# cath-resolve-hits
RUN mkdir -p cath-tools/
# CDD
RUN mkdir -p blast/ncbi-blast-2.15.0+/bin/
RUN mkdir -p rpsbproc/RpsbProc-x64-linux
# Easle
RUN mkdir -p hmmer-3.3/easel/miniapps
# epa-ng
RUN mkdir -p epa-ng/bin/
# HMMER2
RUN mkdir -p hmmer2/bin/
# HMMER3
RUN mkdir -p hmmer3/bin/
# PROSITE pfsearch
RUN mkdir -p pftools/var/lib/pftools/bin/
# PROSITE ps-scan
RUN mkdir -p pftools/var/lib/pftools/bin/

# add scripts that will provide the mock output for the unit tests

# cath-resolve-hits
COPY tests/unit_tests/mock_bins/cath-resolve-hits cath-tools/cath-resolve-hits
# CDD
COPY tests/unit_tests/mock_bins/cdd-blast blast/ncbi-blast-2.15.0+/bin/rpsblast
COPY tests/unit_tests/mock_bins/cdd-rpsbproc rpsbproc/RpsbProc-x64-linux
# Easel
COPY tests/unit_tests/mock_bins/esl-translate hmmer-3.3/easel/miniapps/esl-translate
# epa-ng
COPY tests/unit_tests/mock_bins/epa-ng epa-ng/bin/epa-ng
# HMMER2: hmmpfam
COPY tests/unit_tests/mock_bins/hmmer-2-hmmpfam hmmer2/bin/hmmpfam
# HMMER3: hmmsearch
COPY tests/unit_tests/mock_bins/hmmer-3-hmmsearch hmmer3/bin/hmmsearch
# HMMER3: hmmscan
COPY tests/unit_tests/mock_bins/hmmer-3-hmmscan hmmer3/bin/hmmscan
# PROSITE pf-search
COPY tests/unit_tests/mock_bins/pfsearch pftools/var/lib/pftools/bin/pfsearchV3
# PROSITE ps-scan
COPY tests/unit_tests/mock_bins/ps-scan pftools/var/lib/pftools/bin/ps_scan.pl
# SignalP
WORKDIR /
COPY tests/unit_tests/mock_bins/signalp signalp6

# Copy the files that will form the outputs/inputs to be used in the unit tests
WORKDIR /opt/interproscan6
COPY interproscan/ interproscan/
COPY utilities/ utilities/
COPY tests/ tests/
COPY interproscan.nf interproscan.nf
COPY nextflow.config nextflow.config
COPY README.md README.md
COPY LICENSE LICENSE
