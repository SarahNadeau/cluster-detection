#!/bin/bash
# This script is to run the analysis workflow specified in analyze_boston_confa.nf

set -e

nextflow \
    run \
    analyze_boston_confa.nf \
    -profile docker \
    --context_group_by 'division' \
    --max_context_sequences 100 \
    -resume
