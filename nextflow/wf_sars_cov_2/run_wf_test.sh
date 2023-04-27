#!/bin/bash
# This script is to run test data through the workflow analyze_clusters_sars_cov_2.

set -euo pipefail

# Run the workflow
nextflow \
    run \
    -c ../nextflow.config \
    analyze_clusters_sars_cov_2.nf \
    -profile docker \
    -resume
