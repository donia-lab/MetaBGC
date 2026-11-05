# MetaBGC container: bioinformatics tool-chain + the metabgc Python package.
# Self-contained (does not depend on a prebuilt base image).
FROM ubuntu:20.04
LABEL maintainer="Abhishek Biswas <ab50@princeton.edu>"

ENV DEBIAN_FRONTEND=noninteractive
# Python and Docker are not getting along encoding-wise
ENV LANG=C.UTF-8

# External tools required by metabgc: MUSCLE, HMMER, BLAST+, CD-HIT,
# EMBOSS (transeq) and ART (read simulation).
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
       python3-pip python3-dev \
       g++ rsync zip make git \
       muscle \
       hmmer \
       ncbi-blast+ \
       cd-hit \
       emboss \
       art-nextgen-simulation-tools \
  && ln -sf /usr/bin/python3 /usr/local/bin/python \
  && python3 -m pip install --no-cache-dir --upgrade pip \
  && rm -rf /var/lib/apt/lists/*

# Install metabgc from source.
COPY . /metabgc
RUN python3 -m pip install --no-cache-dir /metabgc

# Thin wrapper that forwards args to the metabgc CLI.
COPY docker-metabgc/run /usr/local/bin/run
RUN chmod +x /usr/local/bin/run

VOLUME ["/input", "/output"]
WORKDIR /output

ENTRYPOINT ["/usr/local/bin/run"]
