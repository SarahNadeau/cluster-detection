#!/bin/bash
# This script is to run the clustertracker workflow specified in clustertracker.nf

set -e

nextflow \
    run \
    clustertracker.nf \
    -profile docker \
    --input_fasta "data/sars_cov_2_boston_confa_outbreak.fasta"
