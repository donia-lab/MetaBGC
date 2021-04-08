# metabgc 2x container
FROM metabgc/release2x:latest
LABEL maintainer="Abhishek Biswas <ab50@princeton.edu>"

# Python and Docker are not getting along encoding-wise
ENV LANG C.UTF-8

# Grab metabgc
COPY . /metabgc

RUN pip3 install /metabgc

ADD docker/run /usr/local/bin/run

VOLUME ["/input", "/output"]
WORKDIR /output

ENTRYPOINT ["/usr/local/bin/run"]
