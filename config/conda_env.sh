#!/bin/bash
# Shared conda activation settings for repo shell entrypoints.

CONDA_ACTIVATE=${CONDA_ACTIVATE:-/scicomp/groups/covlab/SARS2Seq/bin/miniconda/bin/activate}
CONDA_ENV_PREFIX=${CONDA_ENV_PREFIX:-/scicomp/groups/covlab/SARS2Seq/bin/miniconda/envs/prop_model-pure}

activate_prop_model_env() {
    if [[ ! -f "${CONDA_ACTIVATE}" ]]; then
        echo "ERROR: conda activate script not found: ${CONDA_ACTIVATE}" >&2
        echo "Set CONDA_ACTIVATE to a valid path before running this pipeline." >&2
        return 1
    fi

    if [[ ! -d "${CONDA_ENV_PREFIX}" ]]; then
        echo "ERROR: conda environment prefix not found: ${CONDA_ENV_PREFIX}" >&2
        echo "Set CONDA_ENV_PREFIX to a valid env path before running this pipeline." >&2
        return 1
    fi

    source "${CONDA_ACTIVATE}" "${CONDA_ENV_PREFIX}"
}
