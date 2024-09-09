FROM debian:stable-20240904-slim AS downloader

ARG SENTIEON_VERSION
RUN test -n "$SENTIEON_VERSION"

LABEL container.base.image="debian:stable-20240904-slim" \
      software.version="${SENTIEON_VERSION}" \
      software.website="https://www.sentieon.com/"

# Install samtools
RUN apt-get update && apt-get install -y curl bzip2 autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libdeflate-dev && \
    mkdir -p /opt/samtools/ && \
    curl -L "https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2" | \
      tar -C /opt/samtools/ -jxf - && \
    cd /opt/samtools/samtools-1.20 && \
    ./configure && \
    make install

# Install bcftools
RUN apt-get update && apt-get install -y bzip2 autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libperl-dev libgsl0-dev && \
    mkdir -p /opt/bcftools/ && \
    curl -L "https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2" | \
      tar -C /opt/bcftools/ -jxf - && \
    cd /opt/bcftools/bcftools-1.20/ && \
    ./configure && \
    make install

# Install bedtools
RUN apt-get update && apt-get install -y curl && \
    curl -L -o /usr/local/bin/bedtools-2.30.0 "https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary"

# Install igzip
RUN apt-get update && apt-get install -y curl autoconf automake libtool make gcc nasm && \
    mkdir -p /opt/isa-l && \
    curl -L "https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz" | \
        tar -C /opt/isa-l -zxf - && \
    cd /opt/isa-l/isa-l-2.30.0 && \
    ./autogen.sh && \
    ./configure --prefix=/usr --libdir=/usr/lib && \
    make install

# Install mosdepth
RUN apt-get update && apt-get install -y curl && \
    curl -L -o /usr/local/bin/mosdepth-0.3.9 "https://github.com/brentp/mosdepth/releases/download/v0.3.9/mosdepth"

# Download the Sentieon software
RUN apt-get update && apt-get install -y curl && \
    mkdir -p /opt/sentieon/ && \
    curl -L "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${SENTIEON_VERSION}.tar.gz" | \
      tar -zxf - -C /opt/sentieon/

# Install poetry
RUN apt-get update && apt-get install -y git python3 python3-venv && \
    python3 -m venv /opt/poetry-venv && \
    VIRTUAL_ENV=/opt/poetry-venv PATH=/opt/poetry-venv/bin:$PATH pip install poetry

# Build the sentieon-cli
COPY data /opt/sentieon-cli/data
COPY pyproject.toml /opt/sentieon-cli/pyproject.toml
COPY sentieon_cli /opt/sentieon-cli/sentieon_cli
COPY README.md /opt/sentieon-cli/README.md
RUN cd /opt/sentieon-cli/ && \
    /opt/poetry-venv/bin/poetry build -f wheel

# Build the container
FROM debian:stable-20240904-slim
ARG SENTIEON_VERSION
ENV SENTIEON_VERSION=$SENTIEON_VERSION

COPY --from=downloader /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION} /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}
COPY --from=downloader /usr/bin/igzip /usr/bin/igzip
COPY --from=downloader /usr/lib/libisal.a /usr/lib/libisal.a
COPY --from=downloader /usr/lib/libisal.so.2.0.30 /usr/lib/libisal.so.2.0.30
COPY --from=downloader /usr/lib/libisal.la /usr/lib/libisal.la
COPY --from=downloader /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=downloader /usr/local/bin/bcftools /usr/local/bin/bcftools
COPY --from=downloader /usr/local/bin/bedtools-2.30.0 /usr/local/bin/bedtools-2.30.0
COPY --from=downloader /usr/local/bin/mosdepth-0.3.9 /usr/local/bin/mosdepth-0.3.9
COPY --from=downloader /opt/sentieon-cli/dist /opt/sentieon-cli/dist

CMD ["/bin/bash"]

# Create links
RUN cd /usr/local/lib && \
    ln -s libisal.so.2.0.30 libisal.so.2 && \
    ln -s libisal.so.2 libisal.so && \
    cd /usr/local/bin/ && \
    ln -s bedtools-2.30.0 bedtools && \
    ln -s mosdepth-0.3.9 mosdepth && \
    chmod ugo+x bedtools-2.30.0 mosdepth-0.3.9

# Install jemalloc as the recommended memory allocation tool, see https://support.sentieon.com/appnotes/jemalloc/
# Install procps for process monitoring
# Install other dependencies
RUN apt-get update && \
    apt-get install -y libjemalloc2 procps libdeflate-dev libbz2-dev \
      liblzma-dev libcurl4-gnutls-dev libssl-dev libperl-dev libgsl0-dev \
      libncurses5-dev 

ENV SENTIEON_INSTALL_DIR=/opt/sentieon/sentieon-genomics-$SENTIEON_VERSION
ENV PATH $SENTIEON_INSTALL_DIR/bin/:$PATH
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2

# A default jemalloc configuration that should work well for most use-cases, see http://jemalloc.net/jemalloc.3.html
ENV MALLOC_CONF=metadata_thp:auto,background_thread:true,dirty_decay_ms:30000,muzzy_decay_ms:30000

# Create a venv for the sentieon-cli
RUN apt-get update && apt-get install -y git python3 python3-venv curl && \
    python3 -m venv /opt/sentieon-cli-venv
ENV VIRTUAL_ENV /opt/sentieon-cli-venv
ENV PATH /opt/sentieon-cli-venv/bin:$PATH

# Install multiqc into the venv
RUN pip install multiqc

# Install the sentieon-cli into the venv
RUN pip install /opt/sentieon-cli/dist/*.whl

# Test the container
RUN sentieon driver --help && \
    igzip --help && \
    samtools --help && \
    bcftools --help && \
    bedtools --help && \
    multiqc -h && \
    mosdepth -h && \
    sentieon-cli -h
