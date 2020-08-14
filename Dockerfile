FROM ubuntu:20.04

RUN apt update \
  && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends \
  bwa \
  exonerate \
  ncbi-blast+ \
  parallel \
  python3-biopython \
  r-cran-ggplot2 \
  r-cran-reshape2 \
  samtools \
  spades \
  time \
  && rm -rf /var/lib/apt/lists/* \
  && update-alternatives --install /usr/bin/python python /usr/bin/python3 1

COPY *.py gene_recovery_heatmap_ggplot.R /usr/local/bin/

# Avoid conflicting user-installed Python packages when using Singularity
ENV PYTHONNOUSERSITE=1
