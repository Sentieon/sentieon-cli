FROM debian:13.4-slim AS downloader

# Tool versions and pinned SHA-256 checksums. Bumping a version is a
# two-line edit (version + hash) in this block.
ARG SENTIEON_VERSION
ARG SAMTOOLS_VERSION=1.23.1
ARG SAMTOOLS_SHA256=32266198a4bc6a6df395d8526688c9697d9c8e472f888c749fdde2e08ea88dd2
ARG BCFTOOLS_VERSION=1.23.1
ARG BCFTOOLS_SHA256=01899a46f9420cdc1385d52fcfc84cce2806f9c996b787081a90d7dfc85eafa3
ARG BEDTOOLS_VERSION=2.30.0
ARG BEDTOOLS_SHA256=e85d74b6c11b664c05176b1dbf7d2891ad0383ae93805db2d29034db5c2d80ce
ARG ISA_L_VERSION=2.30.0
ARG ISA_L_SHA256=bcf592c04fdfa19e723d2adf53d3e0f4efd5b956bb618fed54a1108d76a6eb56
ARG MOSDEPTH_VERSION=0.3.9
ARG MOSDEPTH_SHA256=ad94aa2bd199c6fd103821d00566765d611f2a4741ce44c5ab723d3b800b911b
ARG VG_VERSION=1.73.0
ARG VG_SHA256=
# Sentieon's KMC fork adds the `-fa /dev/stdin` patch used by the pangenome pipeline.
ARG KMC_VERSION=3.2.4
ARG KMC_RELEASE_TAG=v3.2.4-pipe2
ARG KMC_SHA256=d8bdf8edcd0577dba32e86e8f194b2eb04eb168c75a3e6d46721d5bef515ac96
# Pinned Poetry toolchain; should match the version used to generate poetry.lock.
ARG POETRY_VERSION=2.3.4
ARG POETRY_PLUGIN_EXPORT_VERSION=1.9.0

RUN test -n "$SENTIEON_VERSION"

# Install all build dependencies for the downloader stage in one pass.
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates bzip2 autoconf automake libtool make gcc perl \
      nasm git python3 python3-venv zlib1g-dev libbz2-dev liblzma-dev \
      libcurl4-gnutls-dev libssl-dev libncurses5-dev libdeflate-dev \
      libperl-dev libgsl0-dev && \
    rm -rf /var/lib/apt/lists/*

# Install samtools
RUN mkdir -p /opt/samtools/ && \
    curl -fL -o /tmp/samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    echo "${SAMTOOLS_SHA256}  /tmp/samtools.tar.bz2" | sha256sum -c - && \
    tar -C /opt/samtools/ -jxf /tmp/samtools.tar.bz2 && \
    rm /tmp/samtools.tar.bz2 && \
    cd /opt/samtools/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make install && \
    strip /usr/local/bin/samtools

# Install bcftools
RUN mkdir -p /opt/bcftools/ && \
    curl -fL -o /tmp/bcftools.tar.bz2 "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    echo "${BCFTOOLS_SHA256}  /tmp/bcftools.tar.bz2" | sha256sum -c - && \
    tar -C /opt/bcftools/ -jxf /tmp/bcftools.tar.bz2 && \
    rm /tmp/bcftools.tar.bz2 && \
    cd /opt/bcftools/bcftools-${BCFTOOLS_VERSION}/ && \
    ./configure && \
    make install && \
    strip /usr/local/bin/bcftools

# Install bedtools
RUN curl -fL -o /usr/local/bin/bedtools-${BEDTOOLS_VERSION} "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static.binary" && \
    echo "${BEDTOOLS_SHA256}  /usr/local/bin/bedtools-${BEDTOOLS_VERSION}" | sha256sum -c - && \
    chmod +x /usr/local/bin/bedtools-${BEDTOOLS_VERSION}

# Install igzip
RUN mkdir -p /opt/isa-l && \
    curl -fL -o /tmp/isa-l.tar.gz "https://github.com/intel/isa-l/archive/refs/tags/v${ISA_L_VERSION}.tar.gz" && \
    echo "${ISA_L_SHA256}  /tmp/isa-l.tar.gz" | sha256sum -c - && \
    tar -C /opt/isa-l -zxf /tmp/isa-l.tar.gz && \
    rm /tmp/isa-l.tar.gz && \
    cd /opt/isa-l/isa-l-${ISA_L_VERSION} && \
    ./autogen.sh && \
    ./configure --prefix=/usr --libdir=/usr/lib && \
    make install

# Install mosdepth
RUN curl -fL -o /usr/local/bin/mosdepth-${MOSDEPTH_VERSION} "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth" && \
    echo "${MOSDEPTH_SHA256}  /usr/local/bin/mosdepth-${MOSDEPTH_VERSION}" | sha256sum -c - && \
    chmod +x /usr/local/bin/mosdepth-${MOSDEPTH_VERSION}

# Install vg (variation graph toolkit)
RUN curl -fL -o /usr/local/bin/vg-${VG_VERSION} "https://github.com/vgteam/vg/releases/download/v${VG_VERSION}/vg" && \
    echo "${VG_SHA256}  /usr/local/bin/vg-${VG_VERSION}" | sha256sum -c - && \
    chmod +x /usr/local/bin/vg-${VG_VERSION}

# Install kmc (k-mer counter); Sentieon fork provides the `-fa /dev/stdin` patch
RUN mkdir -p /opt/kmc-${KMC_VERSION} && \
    curl -fL -o /tmp/kmc.tar.gz "https://github.com/Sentieon/KMC/releases/download/${KMC_RELEASE_TAG}/KMC${KMC_VERSION}.linux.x64.tar.gz" && \
    echo "${KMC_SHA256}  /tmp/kmc.tar.gz" | sha256sum -c - && \
    tar -C /opt/kmc-${KMC_VERSION} -zxf /tmp/kmc.tar.gz && \
    rm /tmp/kmc.tar.gz && \
    cp /opt/kmc-${KMC_VERSION}/bin/kmc /usr/local/bin/kmc-${KMC_VERSION} && \
    cp /opt/kmc-${KMC_VERSION}/bin/kmc_tools /usr/local/bin/kmc_tools-${KMC_VERSION} && \
    cp /opt/kmc-${KMC_VERSION}/bin/kmc_dump /usr/local/bin/kmc_dump-${KMC_VERSION} && \
    chmod +x /usr/local/bin/kmc-${KMC_VERSION} \
      /usr/local/bin/kmc_tools-${KMC_VERSION} \
      /usr/local/bin/kmc_dump-${KMC_VERSION}

# Download the Sentieon software
RUN mkdir -p /opt/sentieon/ && \
    curl -fL "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${SENTIEON_VERSION}.tar.gz" | \
      tar -zxf - -C /opt/sentieon/

# Install poetry
RUN python3 -m venv /opt/poetry-venv && \
    /opt/poetry-venv/bin/pip install --no-cache-dir "poetry==${POETRY_VERSION}"

# Build the sentieon-cli
COPY data /opt/sentieon-cli/data
COPY pyproject.toml /opt/sentieon-cli/pyproject.toml
COPY poetry.lock /opt/sentieon-cli/poetry.lock
COPY sentieon_cli /opt/sentieon-cli/sentieon_cli
COPY README.md /opt/sentieon-cli/README.md
RUN /opt/poetry-venv/bin/pip install --no-cache-dir \
      "poetry-plugin-export==${POETRY_PLUGIN_EXPORT_VERSION}" && \
    cd /opt/sentieon-cli/ && \
    /opt/poetry-venv/bin/poetry build -f wheel && \
    /opt/poetry-venv/bin/poetry export \
      --only main,container \
      --format requirements.txt \
      -o /opt/sentieon-cli/dist/requirements.txt

# Build the Python environment in an isolated stage. The venv is
# later COPY'd over as-is.
FROM debian:13.4-slim AS python-builder
ARG SENTIEON_VERSION

COPY --from=downloader /opt/sentieon-cli/dist /opt/sentieon-cli/dist

# Remove the kaleido and pyarrow packages to reduce image size.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      python3 python3-venv python3-dev \
      gcc g++ make zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    python3 -m venv /opt/sentieon-cli-venv && \
    /opt/sentieon-cli-venv/bin/pip install --no-cache-dir \
      -r /opt/sentieon-cli/dist/requirements.txt && \
    /opt/sentieon-cli-venv/bin/pip install --no-cache-dir --no-deps \
      /opt/sentieon-cli/dist/*.whl && \
    /opt/sentieon-cli-venv/bin/pip uninstall -y kaleido pyarrow && \
    /opt/sentieon-cli-venv/bin/pip uninstall -y pip setuptools wheel && \
    find /opt/sentieon-cli-venv -depth -type d -name __pycache__ -exec rm -rf {} +

# Build the container
FROM debian:13.4-slim
ARG SENTIEON_VERSION
ARG BEDTOOLS_VERSION=2.30.0
ARG MOSDEPTH_VERSION=0.3.9
ARG VG_VERSION=1.62.0
ARG KMC_VERSION=3.2.4
ARG VCS_REF
ARG BUILD_DATE

ENV SENTIEON_VERSION=$SENTIEON_VERSION

# OCI image metadata (https://github.com/opencontainers/image-spec/blob/main/annotations.md)
LABEL org.opencontainers.image.title="sentieon-cli" \
      org.opencontainers.image.description="Pipeline implementations for the Sentieon software" \
      org.opencontainers.image.vendor="Sentieon" \
      org.opencontainers.image.url="https://www.sentieon.com/" \
      org.opencontainers.image.source="https://github.com/Sentieon/sentieon-cli" \
      org.opencontainers.image.licenses="BSD-2-Clause" \
      org.opencontainers.image.base.name="debian:13.4-slim" \
      org.opencontainers.image.version="${SENTIEON_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}"

COPY --from=downloader /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION} /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}
COPY --from=downloader /usr/bin/igzip /usr/bin/igzip
COPY --from=downloader /usr/lib/libisal.so.2.0.30 /usr/lib/libisal.so.2.0.30
COPY --from=downloader /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=downloader /usr/local/bin/bcftools /usr/local/bin/bcftools
COPY --from=downloader /usr/local/bin/bedtools-${BEDTOOLS_VERSION} /usr/local/bin/bedtools-${BEDTOOLS_VERSION}
COPY --from=downloader /usr/local/bin/mosdepth-${MOSDEPTH_VERSION} /usr/local/bin/mosdepth-${MOSDEPTH_VERSION}
COPY --from=downloader /usr/local/bin/vg-${VG_VERSION} /usr/local/bin/vg-${VG_VERSION}
COPY --from=downloader /usr/local/bin/kmc-${KMC_VERSION} /usr/local/bin/kmc-${KMC_VERSION}
COPY --from=downloader /usr/local/bin/kmc_tools-${KMC_VERSION} /usr/local/bin/kmc_tools-${KMC_VERSION}
COPY --from=downloader /usr/local/bin/kmc_dump-${KMC_VERSION} /usr/local/bin/kmc_dump-${KMC_VERSION}
COPY --from=downloader /opt/sentieon-cli/dist /opt/sentieon-cli/dist

# Create symlinks for libisal (the .so lives on the default linker path in
# /usr/lib) and for the unversioned names of the third-party tools.
RUN cd /usr/lib && \
    ln -s libisal.so.2.0.30 libisal.so.2 && \
    ln -s libisal.so.2 libisal.so && \
    cd /usr/local/bin/ && \
    ln -s bedtools-${BEDTOOLS_VERSION} bedtools && \
    ln -s mosdepth-${MOSDEPTH_VERSION} mosdepth && \
    ln -s vg-${VG_VERSION} vg && \
    ln -s kmc-${KMC_VERSION} kmc && \
    ln -s kmc_tools-${KMC_VERSION} kmc_tools && \
    ln -s kmc_dump-${KMC_VERSION} kmc_dump

# Install runtime dependencies.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      libjemalloc2 procps libdeflate0 libbz2-1.0 liblzma5 \
      libcurl3-gnutls libssl3 libperl5.40 libgsl28 libncurses6 \
      curl python3 && \
    rm -rf /var/lib/apt/lists/*

ENV SENTIEON_INSTALL_DIR=/opt/sentieon/sentieon-genomics-$SENTIEON_VERSION
ENV PATH=$SENTIEON_INSTALL_DIR/bin/:$PATH
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2

# A default jemalloc configuration that should work well for most use-cases, see http://jemalloc.net/jemalloc.3.html
ENV MALLOC_CONF=metadata_thp:auto,background_thread:true,dirty_decay_ms:30000,muzzy_decay_ms:30000

# Copy the pre-built venv from the python-builder stage. No compilers
# are installed in this stage.
COPY --from=python-builder /opt/sentieon-cli-venv /opt/sentieon-cli-venv
ENV VIRTUAL_ENV=/opt/sentieon-cli-venv
ENV PATH=/opt/sentieon-cli-venv/bin:$PATH

# Create a non-root user for running the pipelines
RUN useradd --create-home --uid 1001 --shell /bin/bash sentieon
USER sentieon
WORKDIR /home/sentieon

# Test the container
RUN sentieon driver --help && \
    igzip --help && \
    samtools --help && \
    bcftools --help && \
    bedtools --help && \
    multiqc -h && \
    mosdepth -h && \
    vg version && \
    kmc --help && \
    sentieon-cli -h
