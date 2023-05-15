#!/bin/bash
# This script is to run the analysis workflow specified in analyze_boston_confa.nf on the gizmo cluster.

set -euo pipefail

# Load modules (specific to rhino/gizmo)
ml purge
ml nextflow/23.04.0
ml Apptainer/1.1.6

# Up the recursion limit for augur for big tree
export AUGUR_RECURSION_LIMIT=10000

# Run the workflow
nextflow \
    run \
    -c ../nextflow.gizmo.notower.config \
    -N snadeau@fredhutch.org \
    analyze_clusters_sars_cov_2.nf \
    -profile apptainer \
    --input_fasta "../../clean_data/sars_cov_2_lemiux_boston/outbreak_samples_confa.fasta" \
    --input_metadata "../../clean_data/sars_cov_2_lemiux_boston/clustertracker_metadata.txt" \
    --input_metadata_nextstrain "../../clean_data/sars_cov_2_lemiux_boston/nextstrain_metadata.txt" \
    --trait_name "division" \
    --context_region_name '.' \
    --reference_fasta "../../clean_data/sars_cov_2_lemiux_boston/reference_NC_045512v2.fa" \
    --filter_similarity_specs "--min-length 20000 --max-date '2020-07-01' --subsample-max-sequences 100" \
    --filter_geocontext_specs "--min-length 20000 --max-date '2020-07-01' --group-by region year month --subsample-max-sequences 200" \
    --outgroup_taxon "NC_045512v2" \
    --nextstrain_refine_params "--coalescent opt --date-confidence --clock-rate 0.0008 --clock-std-dev 0.004 --root NC_045512v2" \
    -with-report sars_cov_2_lemiux_boston_benchmark.html \
    -stub-run
