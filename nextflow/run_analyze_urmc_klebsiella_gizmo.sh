#!/bin/bash
# This script is to run the analysis workflow specified in analyze_urmc_klebsiella.nf on the gizmo cluster.

set -euo pipefail

# Load modules (specific to rhino/gizmo)
ml purge
ml nextflow/23.04.0
ml Apptainer/1.1.6

# Run the workflow
nextflow \
    run \
    -c nextflow.gizmo.config \
    analyze_urmc_klebsiella.nf \
    -profile apptainer \
    -N snadeau@fredhutch.org \
    --input_fasta_dir '../clean_data/all_samples_for_analysis' \
    -resume
