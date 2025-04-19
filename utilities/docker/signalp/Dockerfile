FROM python:3.7-slim
ENV TZ=Europe/London
RUN apt-get update -y && apt-get install -y procps
RUN pip install "matplotlib>3.3.2" "numpy>1.19.2" "torch>1.7.0,<2" "tqdm>4.46.1"
