FROM python:3.9-slim AS builder

ENV TZ=Europe/London

WORKDIR /opt/tmbed
COPY . /opt/tmbed

RUN pip install --no-cache-dir numpy==2.0.2 \
 && pip install --no-cache-dir h5py==3.13.0 \
 && pip install --no-cache-dir sentencepiece==0.2.0 \
 && pip install --no-cache-dir tqdm==4.67.1 \
 && pip install --no-cache-dir transformers==4.51.3 \
 && pip install --no-cache-dir typer==0.15.3 \
 && pip install --no-cache-dir torch==2.7.0+cu128 --extra-index-url https://download.pytorch.org/whl/cu128


# Download and install TMbed v1.0.2
ADD https://github.com/BernhoferM/TMbed/archive/refs/tags/v1.0.2.tar.gz /opt/tmbed/tmbed.tar.gz
RUN tar -xzf tmbed.tar.gz \
 && cd TMbed-1.0.2 \
 && pip install --no-cache-dir .

# Download the models
RUN tmbed download -m /opt/tmbed/models

FROM python:3.9-slim

ENV TZ=Europe/London

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

# Copy Python packages and bins
COPY --from=builder /usr/local/lib/python3.9 /usr/local/lib/python3.9
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy the downloaded models
COPY --from=builder /opt/tmbed/models /opt/tmbed/models

ENV CUBLAS_WORKSPACE_CONFIG=:4096:8
