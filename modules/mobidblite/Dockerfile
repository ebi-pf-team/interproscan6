FROM ubuntu:24.04

LABEL name="interpro/idrpred"
LABEL tag="1.0.3"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/venv/bin:/root/.local/bin:${PATH}"
ENV UV_PYTHON_INSTALL_DIR="/opt/uv/python"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        procps \
        unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    uv python install 3.12 && \
    uv venv --python 3.12 /opt/venv && \
    curl -L -o idrpred.zip https://github.com/matthiasblum/idrpred/archive/refs/tags/1.0.3.zip && \
    unzip idrpred.zip && \
    uv pip install /opt/idrpred-1.0.3 && \
    rm -rf /opt/idrpred-1.0.3 /opt/idrpred.zip
