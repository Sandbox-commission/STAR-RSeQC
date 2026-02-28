# STAR-RSeQC: Complete Bioinformatics Pipeline in a Container
# Multi-stage build: compile in builder, minimal runtime image

FROM ubuntu:22.04 as builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    build-essential \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Copy source code
WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src/ ./src/

# Build release binary
RUN cargo build --release

# ─────────────────────────────────────────────────────────────────────────────

FROM mambaorg/micromamba:latest

# Set working directory
WORKDIR /data

# Install system dependencies
RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment for STAR
ENV MAMBA_DOCKERFILE_WORKDIR=/data
RUN micromamba create -y -n star \
    -c bioconda \
    -c conda-forge \
    star=2.7.11b \
    samtools \
    && micromamba clean --all -y

# Create conda environment for RSeQC
RUN micromamba create -y -n rseqc \
    -c bioconda \
    -c conda-forge \
    rseqc \
    python \
    && micromamba clean --all -y

# Copy pre-built binary from builder
COPY --from=builder /build/target/release/star-rseqc /usr/local/bin/star-rseqc

# Create config directory
RUN mkdir -p /root/.config/star-rseqc

# Set up environment activation
RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'set -e' >> /entrypoint.sh && \
    echo 'eval "$(micromamba shell hook --shell bash)"' >> /entrypoint.sh && \
    echo 'export PATH="/opt/conda/envs/star/bin:/opt/conda/envs/rseqc/bin:$PATH"' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

# Default command
ENTRYPOINT ["/entrypoint.sh"]
CMD ["star-rseqc", "-h"]

# Labels
LABEL maintainer="STAR-RSeQC"
LABEL description="STAR 2-pass alignment + RSeQC QC pipeline"
LABEL version="0.1.0"
