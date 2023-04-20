#!/bin/bash
# This script is to run the nextstrain workflow specified in nextstrain_filter.nf

set -e

nextflow \
    run \
    nextstrain_filter.nf \
    -profile docker \
    -resume \
    --region_name "north-america"
