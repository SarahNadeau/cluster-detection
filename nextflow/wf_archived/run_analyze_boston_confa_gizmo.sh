#!/bin/bash
# This script is to run the analysis workflow specified in analyze_boston_confa.nf on the gizmo cluster.

set -euo pipefail

# Load modules (specific to rhino/gizmo)
ml purge
ml nextflow/23.04.0
ml Apptainer/1.1.6

# Run the workflow
nextflow \
    run \
    -c nextflow.gizmo.config \
    analyze_boston_confa.nf \
    -profile apptainer \
    --context_group_by 'region' \
    --max_context_sequences 100 \
    -resume
