# SC2 Proportion Modeling - Multi-Stage Docker Build
# Stage 1: Base Python environment
# Stage 2: R dependencies via Conda
# Stage 3: Python dependencies via Poetry
# Stage 4: Runtime optimized image

# ============================================================================
# Stage 1: Base - Python 3.11 + System Dependencies
# ============================================================================
FROM python:3.11-slim as base

LABEL maintainer="CDC Emerging Infectious Diseases Branch"
LABEL version="2.0.0"
LABEL description="SC2 Proportion Modeling Pipeline - Docker Container"

# Set non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    wget \
    git \
    ca-certificates \
    libssl-dev \
    libffi-dev \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip, setuptools, wheel
RUN python -m pip install --upgrade pip setuptools wheel

# ============================================================================
# Stage 2: R Dependencies via Conda
# ============================================================================
FROM base as r-deps

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean -afy

# Set PATH to include conda
ENV PATH="/opt/conda/bin:${PATH}"

# Add conda-forge and create R environment
RUN /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda create -n r-env -y \
    r-base \
    r-survey \
    r-ggplot2 \
    r-dplyr \
    r-tidyr \
    r-gridextra \
    r-scales && \
    /opt/conda/bin/conda clean -afy

ENV R_HOME=/opt/conda/envs/r-env \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# ============================================================================
# Stage 3: Python Dependencies via Poetry
# ============================================================================
FROM r-deps as builder

# Install Poetry
RUN python -m pip install poetry>=1.5.0

# Set working directory
WORKDIR /build

# Copy poetry files
COPY pyproject.toml poetry.lock* ./

# Install Python dependencies
RUN poetry config virtualenvs.create false && \
    poetry install --no-dev --no-interaction --no-ansi

# ============================================================================
# Stage 4: Runtime - Optimized Final Image
# ============================================================================
FROM r-deps as runtime

# Set non-interactive mode (inherit from base)
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PATH="/opt/conda/bin:${PATH}" \
    R_HOME=/opt/conda/envs/r-env \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# Install runtime Python packages only
RUN python -m pip install --no-cache-dir \
    polars>=0.20.0 \
    pydantic>=2.0.0 \
    pydantic-settings>=2.0.0 \
    pyyaml>=6.0 \
    sqlalchemy>=2.0.0 \
    pandas>=2.0.0 \
    numpy>=1.24.0 \
    cmdstanpy>=1.1.0 \
    rpy2>=3.5.0 \
    pyarrow>=12.0.0 \
    requests>=2.31.0 \
    python-dotenv>=1.0.0 \
    loguru>=0.7.0

# Create application directories
RUN mkdir -p /app /app/data /app/config /app/results /app/cache /app/logs

# Create non-root user for security
RUN groupadd -r appuser && \
    useradd -r -g appuser -u 1000 appuser

# Copy application code from builder stage
COPY --chown=appuser:appuser src/ /app/src/
COPY --chown=appuser:appuser tests/ /app/tests/
COPY --chown=appuser:appuser pyproject.toml /app/

# Install application in editable mode
WORKDIR /app
RUN python -m pip install -e . --no-deps

# Fix permissions
RUN chown -R appuser:appuser /app

# Switch to non-root user
USER appuser

# Healthcheck
HEALTHCHECK --interval=60s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import sc2; print('OK')" || exit 1

# Set default command
ENTRYPOINT ["sc2-run"]
CMD ["--help"]

# Labels for tracking
LABEL org.opencontainers.image.version="2.0.0"
LABEL org.opencontainers.image.created="2026-03-18"
LABEL org.opencontainers.image.source="https://github.com/anthropics/sc2-proportion-modeling"
