#!/usr/bin/env bash
# Create or update conda env with all pipeline tools
# Usage: ./bin/setup_env.sh [env_name]
# Example: ./bin/setup_env.sh hw3

set -euo pipefail

ENV_NAME="${1:-hw3}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

if ! command -v conda &>/dev/null; then
    echo "Error: conda not found. Install miniconda first."
    exit 1
fi

echo "Setting up conda env: ${ENV_NAME}"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Updating existing env..."
    conda env update -n "${ENV_NAME}" -f "${ROOT}/envs/hw3_all.yml" --prune
else
    echo "Creating new env..."
    conda env create -n "${ENV_NAME}" -f "${ROOT}/envs/hw3_all.yml"
fi

echo ""
echo "Done. Run pipeline with:"
echo "  nextflow run main.nf -profile local --conda_env \$(conda info --base)/envs/${ENV_NAME} \\"
echo "    --reads '/path/to/*_{1,2}.fastq' \\"
echo "    --reference /path/to/ref.fasta"
