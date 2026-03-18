# syntax=docker/dockerfile:1.7

FROM rust:1.76-bookworm AS builder
WORKDIR /app
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release

FROM mambaorg/micromamba:2.3.3
USER root

# Install bioinformatics runtime dependencies in a dedicated conda env.
RUN micromamba create -y -n bio -c conda-forge -c bioconda \
    star=2.7.11b \
    rseqc=5.0.1 \
    deeptools=3.5.6 \
    samtools=1.20 \
    r-base=4.3 \
    r-ggplot2 && \
    micromamba clean --all --yes

COPY --from=builder /app/target/release/star-rseqc /usr/local/bin/star-rseqc
RUN chmod +x /usr/local/bin/star-rseqc

ENV STAR_RSEQC_STAR_ENV=/opt/conda/envs/bio \
    STAR_RSEQC_RSEQC_ENV=/opt/conda/envs/bio \
    STAR_RSEQC_DEEPTOOLS_ENV=/opt/conda/envs/bio \
    PATH=/opt/conda/envs/bio/bin:${PATH}

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/star-rseqc"]
CMD ["--help"]
