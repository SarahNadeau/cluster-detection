#!/bin/bash

set -euo pipefail

# Load modules (specific to rhino/gizmo)
ml purge
ml nextflow/23.04.0
ml Apptainer/1.1.6

# Run the workflow
nextflow \
    run \
    -c ../nextflow.gizmo.config \
    -N snadeau@fredhutch.org \
    analyze_urmc_klebsiella.nf \
    -profile apptainer \
    -resume
