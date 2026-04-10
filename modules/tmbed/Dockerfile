FROM python:3.9-slim AS builder

LABEL name="interpro/tmbed"
LABEL tag="1.0.2"

ENV TZ=Europe/London

WORKDIR /opt/tmbed

RUN apt-get update && \
    apt-get install -y curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    numpy==2.0.2 \
    h5py==3.13.0 \
    sentencepiece==0.2.0 \
    tqdm==4.67.1 \
    transformers==4.51.3 \
    typer==0.15.3 \
    torch==2.7.0+cu128 --extra-index-url https://download.pytorch.org/whl/cu128

RUN curl -L -O https://github.com/BernhoferM/TMbed/archive/refs/tags/v1.0.2.tar.gz && \
    tar -xzf v1.0.2.tar.gz && \
    cd TMbed-1.0.2 && \
    pip install --no-cache-dir .

# Download the models
RUN tmbed download -m /opt/tmbed/models

FROM python:3.9-slim

ENV TZ=Europe/London

RUN apt-get update && \
    apt-get install -y procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy Python packages and bins
COPY --from=builder /usr/local/lib/python3.9 /usr/local/lib/python3.9
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy the downloaded models
COPY --from=builder /opt/tmbed/models /opt/tmbed/models

ENV CUBLAS_WORKSPACE_CONFIG=:4096:8
