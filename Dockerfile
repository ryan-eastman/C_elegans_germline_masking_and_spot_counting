# GPU container for the RTX 5090 (Blackwell sm_120) — build on the workstation,
# then convert to Apptainer for the HPC (see apptainer.def).
#   docker build -t germquant:gpu .
FROM nvidia/cuda:12.8.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.11 python3.11-venv python3-pip git libgl1 libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*
RUN ln -sf /usr/bin/python3.11 /usr/bin/python && python -m pip install --upgrade pip

WORKDIR /opt/germquant
COPY pyproject.toml README.md ./
COPY src ./src

# Explicit cu128 torch wheel FIRST (do not let resolvers pick a CPU/older build),
# then the package + GPU/SC extras.
RUN pip install --index-url https://download.pytorch.org/whl/cu128 "torch>=2.7" \
    && pip install ".[gpu,sc]"

# fail-fast: the lock is only valid if the GPU is actually sm_120
RUN python -c "import germquant, skan, cellpose; print('germquant ok')"

ENTRYPOINT ["germquant"]
CMD ["--help"]
