#!/bin/bash

set -euo pipefail

# Load modules (specific to rhino/gizmo)
ml purge
ml nextflow/23.04.0
ml Apptainer/1.1.6

# Run the workflow
nextflow \
    run \
    -c ../nextflow.gizmo.notower.config \
    -N snadeau@fredhutch.org \
    analyze_clusters_bacteria.nf \
    -profile apptainer \
    --input_fasta_dir "../../clean_data/klebsiella_malek_urmc/all_samples_for_analysis" \
    --input_metadata "../../clean_data/klebsiella_malek_urmc/clustertracker_metadata.txt" \
    --input_metadata_nextstrain "../../clean_data/klebsiella_malek_urmc/nextstrain_metadata.csv" \
    --reference_fasta "../../clean_data/klebsiella_malek_urmc/reference_NC_015663v1.fna" \
    --trait_name "location" \
    --nextstrain_refine_params "--coalescent opt --root reference_NC_015663v1.fna.ref" \
    -resume \
    -stub-run  # using cached results from 28.04.2023
