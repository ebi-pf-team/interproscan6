# docker build -t interpro/groovy:4.0.27-1 -f docker/groovy.Dockerfile .

FROM groovy:4.0.27-jdk17 AS groovy-runtime
FROM ubuntu:24.04

WORKDIR /opt/interproscan6

COPY --from=groovy-runtime /opt/groovy /opt/groovy
COPY --from=groovy-runtime /opt/java/openjdk /opt/java/openjdk
ENV GROOVY_HOME=/opt/groovy
ENV JAVA_HOME=/opt/java/openjdk
ENV PATH="/opt/groovy/bin:/opt/java/openjdk/bin:$PATH"

COPY lib /opt/interproscan6/lib

RUN apt-get update && \
    apt-get install -y procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
