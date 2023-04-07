#!/bin/bash
# This script is to run the clustertracker workflow specified in clustertracker.nf on the gizmo cluster.

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
    -C ${NXF_CONFIG} \
    run \
    clustertracker.nf \
    -with-report nextflow.report.html \
    -profile singularity

