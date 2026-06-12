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
# config/ and workflow/ are needed at runtime: `germquant run --config config/config.yaml`
# loads the channel map by a path relative to the config dir, and the Snakefile lives here.
COPY config ./config
COPY workflow ./workflow

# Explicit cu128 torch wheel FIRST (do not let resolvers pick a CPU/older build),
# then the package + GPU/SC/workflow extras.
RUN pip install --index-url https://download.pytorch.org/whl/cu128 "torch>=2.7" \
    && pip install ".[gpu,sc,workflow]"

# build-time smoke check (imports only — a GPU is not attached during `docker build`).
RUN python -c "import germquant, skan, cellpose; print('germquant ok')"

# Runtime sm_120 assertion lives in the entrypoint, not the build: `docker build` has no
# GPU, so the (12,0) capability check can only run with `--gpus all` at `docker run`.
# Call `germquant check-gpu` (or the pixi `check-gpu` task) on the 5090/HPC before trusting a run.

ENTRYPOINT ["germquant"]
CMD ["--help"]
