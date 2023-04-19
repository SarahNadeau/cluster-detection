#!/bin/bash
# This script is to run the analysis workflow specified in analyze_boston_confa.nf on the gizmo cluster.

set -euo pipefail

# Load modules (specific to rhino/gizmo)
ml purge
ml nextflow/23.04.0
ml Apptainer/1.1.6

# Run the workflow
echo "WARNING: stub-run is on!"
nextflow \
    run \
    -c nextflow.gizmo.config \
    -N snadeau@fredhutch.org \
    analyze_boston_confa.nf \
    -profile apptainer \
    --context_region_name '.' \
    --context_group_by 'region' \
    --max_context_sequences 500 \
    -resume \
    -stub-run
