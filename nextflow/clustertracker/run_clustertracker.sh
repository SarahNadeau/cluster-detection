#!/bin/bash
# This script is to run the clustertracker workflow specified in clustertracker.nf

set -e

nextflow \
    run \
    clustertracker.nf \
    -profile docker \
    --pb "/Users/snadeau/Documents/repos/cluster-detection/nextflow/clustertracker/assets/test_data/public-2021-06-09.all.masked.nextclade.pangolin.pb.gz"
