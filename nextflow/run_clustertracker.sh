#!/bin/bash
# This script is to run the clustertracker workflow specified in clustertracker.nf

set -e

nextflow \
    run \
    clustertracker.nf \
    -profile docker \
    -resume
