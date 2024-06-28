FROM python:3.7-slim
ENV TZ=Europe/London
WORKDIR /opt/signalp6
COPY . .
RUN apt-get update -y && apt-get upgrade -y && apt-get install -y procps
RUN pip install --no-cache-dir . && \
    SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" ) && \
    cp -r models/* $SIGNALP_DIR/model_weights/
