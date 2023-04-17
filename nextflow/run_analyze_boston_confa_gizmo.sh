#!/bin/bash
# This script is to run the analysis workflow specified in analyze_boston_confa.nf on the gizmo cluster.

set -euo pipefail

# Nextflow Configuration File
NXF_CONFIG=./nextflow.gizmo.config

# Specify a TMPDIR (otherwise apptainer complains it can't write)
export TMPDIR=/fh/scratch/delete30/matsen_e/snadeau/tmp

# Load modules (specific to rhino/gizmo)
ml purge
ml nextflow/23.04.0
ml Apptainer/1.1.6

# Run the workflow
nextflow \
    run \
    -c $NXF_CONFIG \
    analyze_boston_confa.nf \
    -profile apptainer \
    --context_group_by 'region' \
    --max_context_sequences 100 \
    -resume
