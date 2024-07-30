FROM python:3.7-slim
ENV TZ=Europe/London
WORKDIR /opt/phobius
COPY . .
RUN apt-get update -y && apt-get upgrade -y && apt-get install -y procps perl && \
    sed -i "s|/decodeanhmm|/decodeanhmm.64bit|" /opt/phobius/phobius.pl && \
    ln -s /opt/phobius/phobius.pl /usr/local/bin/phobius
