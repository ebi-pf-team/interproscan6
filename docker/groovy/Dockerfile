# Build from the root of the repo (replace 'x')
# docker build -t interpro/groovy:4.0.27-x -f docker/groovy/Dockerfile .

FROM groovy:4.0.27-jdk17 AS groovy-runtime
FROM ubuntu:24.04

WORKDIR /opt/interproscan6

COPY --from=groovy-runtime /opt/groovy /opt/groovy
COPY --from=groovy-runtime /opt/java/openjdk /opt/java/openjdk
ENV GROOVY_HOME=/opt/groovy
ENV JAVA_HOME=/opt/java/openjdk
ENV PATH="/opt/groovy/bin:/opt/java/openjdk/bin:$PATH"

COPY lib /opt/interproscan6/lib
COPY docker/groovy/*groovy /opt/interproscan6

RUN apt-get update && \
    apt-get install -y procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
