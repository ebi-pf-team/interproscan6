FROM ubuntu:latest

RUN apt-get update && \
apt-get upgrade -y && \
DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC \
apt-get install -y autoconf automake autotools-dev bison build-essential cmake flex libtool zlib1g-dev git python3 python3-pip python3-venv

# Create a Python virtual environment
RUN python3 -m venv /opt/venv
# Activate the virtual environment
ENV PATH="/opt/venv/bin:$PATH"

# Install epa-ng
RUN git clone https://github.com/pierrebarbera/epa-ng
RUN cd epa-ng && make 
# the resulting epa-ng binary is located at bin/epa-ng

# Install biopython
RUN pip install biopython
