# ---- stage 1: build PoissonRecon ----
FROM ubuntu:22.04 AS poisson-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    make \
    g++ \
    libjpeg-dev \
    libpng-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/mkazhdan/PoissonRecon.git /build/PoissonRecon

RUN make -j1 -C /build/PoissonRecon poissonrecon COMPILER=gcc


# ---- stage 2: runtime ----
FROM ubuntu:22.04

ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-dev \
    python3-pip \
    libjpeg-dev \
    libpng-dev \
    zlib1g-dev \
    libgomp1 \
    libgl1 \
    libglib2.0-0 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=poisson-builder /build/PoissonRecon/Bin/Linux/PoissonRecon /usr/local/bin/PoissonRecon
RUN chmod +x /usr/local/bin/PoissonRecon

WORKDIR /npet2

COPY pyproject.toml .
COPY libnpet/ libnpet/

RUN pip3 install --no-cache-dir --upgrade pip setuptools && \
    pip3 install --no-cache-dir .

ENV NPET2_ROOT=/data
ENV NPET2_RUNS_ROOT=/data/runs
ENV NPET2_CACHE_ROOT=/data/cache
ENV NPET2_POISSON_RECON_BIN=/usr/local/bin/PoissonRecon
ENV NPET2_RIBOXYZ_API_URL=https://api.ribosome.xyz

VOLUME ["/data"]

ENTRYPOINT ["npet2"]
CMD ["--help"]