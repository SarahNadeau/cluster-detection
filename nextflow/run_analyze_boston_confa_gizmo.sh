#!/bin/bash
# This script is to run the analysis workflow specified in analyze_boston_confa.nf on the gizmo cluster.

set -euo pipefail

# Nextflow Configuration File
NXF_CONFIG=~/nextflow/nextflow.gizmo.config

# Load the Nextflow module (if running on rhino/gizmo)
ml Nextflow/22.10.6

# Load the Singularity module (if running on rhino/gizmo with Singularity)
ml Singularity
# Make sure that the singularity executables are in the PATH
export PATH=$SINGULARITYROOT/bin/:$PATH

# Run the workflow
nextflow \
    run \
    analyze_boston_confa.nf \
    -profile singularity \
    --context_region_name \
    --context_group_by 'region' \
    --max_context_sequences 100 \
    -resume
