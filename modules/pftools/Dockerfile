FROM sibswiss/pftools AS pftools
FROM ubuntu:24.04

LABEL name="interpro/pftools"
LABEL tag="3.2.12"

ENV PATH="/opt/pftools:${PATH}"

COPY --from=pftools /var/lib/pftools/bin/pfsearchV3 /opt/pftools/pfsearchV3
COPY --from=pftools /var/lib/pftools/bin/pfscanV3 /opt/pftools/pfscanV3

RUN apt-get update && \
    apt-get install -y libdata-dump-perl procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*