#!/usr/bin/env bash
# Build and push all HW3 Docker images to Docker Hub
# Usage: ./build_and_push.sh [docker_hub_username]

set -euo pipefail

REGISTRY="${1:-eugeniaarth}"
TAG="1.0"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "Building images for ${REGISTRY}..."

docker build -t "${REGISTRY}/hw3-fastqc:${TAG}" "${ROOT}/docker/fastqc"
docker build -t "${REGISTRY}/hw3-fastp:${TAG}" "${ROOT}/docker/fastp"
docker build -t "${REGISTRY}/hw3-bwa-samtools:${TAG}" "${ROOT}/docker/bwa_samtools"
docker build -f "${ROOT}/docker/python_plot/Dockerfile" -t "${REGISTRY}/hw3-python-plot:${TAG}" "${ROOT}"
docker build -t "${REGISTRY}/hw3-bcftools:${TAG}" "${ROOT}/docker/bcftools"

echo "Pushing images..."
docker push "${REGISTRY}/hw3-fastqc:${TAG}"
docker push "${REGISTRY}/hw3-fastp:${TAG}"
docker push "${REGISTRY}/hw3-bwa-samtools:${TAG}"
docker push "${REGISTRY}/hw3-python-plot:${TAG}"
docker push "${REGISTRY}/hw3-bcftools:${TAG}"

echo "Done."
