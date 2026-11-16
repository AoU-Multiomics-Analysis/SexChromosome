# syntax=docker/dockerfile:1.7
FROM mambaorg/micromamba:2.10.0

LABEL org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/SexChromosome" \
      org.opencontainers.image.title="XYalign for Terra" \
      org.opencontainers.image.version="1.1.6" \
      org.opencontainers.image.description="Unmodified XYalign 1.1.6 with a typed-command wrapper for Terra"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml --strict-channel-priority \
    && micromamba clean --all --yes

# Keep upstream source unchanged. Its Python 2 runtime is supplied by conda.
USER root
ADD --checksum=sha256:f011e8ea3f8771741b707240a02a37dc3829e05b8663f407c8face108614ed80 \
    https://codeload.github.com/SexChrLab/XYalign/tar.gz/a68f484cd256dec37abb1166bb3a8a0fb994351f \
    /tmp/xyalign.tar.gz
ENV PATH=/opt/conda/bin:$PATH \
    MPLBACKEND=Agg \
    PYTHONUNBUFFERED=1
RUN python -m pip install --no-deps --no-build-isolation --no-cache-dir /tmp/xyalign.tar.gz \
    && python -c "from xyalign._version import __version__; assert __version__ == '1.1.6', __version__" \
    && rm /tmp/xyalign.tar.gz /tmp/environment.yml

COPY scripts/xyalign_runner.py /opt/sexchromosome/scripts/xyalign_runner.py
WORKDIR /work
CMD ["python", "/opt/sexchromosome/scripts/xyalign_runner.py", "--help"]
