#!/bin/bash
# This script is to run test data through the workflow analyze_clusters_bacteria.nf

set -euo pipefail

# Run the workflow
nextflow \
    run \
    -c ../nextflow.config \
    analyze_clusters_bacteria.nf \
    -profile docker \
    -with-report analyze_clusters_bacteria.html
