FROM python:3.8-slim
ENV TZ=Europe/London
WORKDIR /opt/mobidb

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y procps git
RUN pip install git+https://github.com/matthiasblum/idrpred
