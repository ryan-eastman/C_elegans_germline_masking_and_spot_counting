# GPU container — the cu128 torch wheel covers the workstation RTX 5090 (Blackwell sm_120) AND the
# HPC A100 (sm_80) / L40 (sm_89), so this single image is PORTABLE across the desktop and RMACC
# Alpine. Build on the workstation, convert to Apptainer for the HPC (see apptainer.def, docs/HPC_ALPINE.md).
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
COPY constraints ./constraints
# the built-in system packs (docs/PACKS.md): found under $GERMQUANT_CONFIG_ROOT/packs, set below
COPY packs ./packs

# Pinned to the validated workstation environment (constraints/desktop-2026-09.txt) so a container
# run resolves the same package versions the goldens were made with. Explicit cu128 torch wheel FIRST
# (do not let resolvers pick a CPU/older build), then the package with the GPU / spots / SC-tracing /
# workflow extras.
RUN pip install --index-url https://download.pytorch.org/whl/cu128 "torch==2.11.0+cu128" "torchvision==0.26.0" \
    && pip install -c constraints/desktop-2026-09.txt ".[gpu,sc,workflow]" spotmax cellacdc cupy-cuda12x

# profiles are found from this checkout when the working directory is elsewhere (`--profile NAME`)
ENV GERMQUANT_CONFIG_ROOT=/opt/germquant

# build-time smoke check (imports only — a GPU is not attached during `docker build`).
RUN python -c "import germquant, cellpose; print('germquant ok')"
# the built-in packs are listed (gate G13); the engine is installed non-editable here, so they are found
# through GERMQUANT_CONFIG_ROOT, and this torch-free entry fails the build if none is found
RUN python -m germquant.packs list && python -c "import germquant.packs as P, sys; sys.exit(0 if P.list_packs() else 1)"

# Runtime sm_120 assertion lives in the entrypoint, not the build: `docker build` has no
# GPU, so the (12,0) capability check can only run with `--gpus all` at `docker run`.
# Call `germquant check-gpu` (or the pixi `check-gpu` task) on the 5090/HPC before trusting a run.

ENTRYPOINT ["germquant"]
CMD ["--help"]
