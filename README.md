# Molecular transmission cluster detection project

This project aims to be a methods guide and resource for molecular transmission cluster detection for public health.
<!-- The full write-up is available at TODO. -->

## Methods implemented:
* HIV-TRACE
* ClusterTracker (matUtils introduce)
* Nextstrain's _augur_

## Resources:
* [clean_data](./clean_data/): Collated metadata for one bacterial, one viral outbreak with ground-truth clustering
* [modules](./nextflow/modules/): Re-usable nextflow processes implementing the different cluster detection methods 
* [wf_bacteria](./nextflow/wf_bacteria/): A nextflow workflow to run the three methods on the bacterial outbreak data set
* [wf_sars_cov_2](./nextflow/wf_sars_cov_2/): A nextflow workflow to run the three methods on the SARS-CoV-2 outbreak data set

## Running workflows

Here we'll run some test data through the SARS-CoV-2 and bacterial cluster detection workflows. 
First, install [Docker](https://docs.docker.com/get-docker/) and [Nextflow](https://www.nextflow.io/).

Then run each workflow:
```
# SARS-CoV-2 test data
cd nextflow/wf_sars_cov_2
nextflow \
    run \
    -c ../nextflow.config \
    analyze_clusters_sars_cov_2.nf \
    -profile docker

# PhiX bacterial test data
cd ../wf_bacteria
nextflow \
    run \
    -c ../nextflow.config \
    analyze_clusters_bacteria.nf \
    -profile docker
```

## Ananlyzing your own data

Each workflow requires the following inputs:
* A reference sequence for alignment and to root the phylogeny
* Focal outbreak genome sequences
    * SARS-CoV-2: FASTA-format assemblies, all stored in the same file
    * Bacteria: FASTA-format assemblies stored in a directory, 1 assembly per file
* Context genome sequences for comparison, to represent generally circulating strains
    * SARS-CoV-2: downloaded by the workflow using Nextstrain's data resources
    * Bacteria: FASTA-format assemblies stored in the same directory as outbreak sequences, 1 assembly per file
* Metadata
    * Two metadata files required, one in ClusterTracker format, one in Nextstrain format, see below for examples
    * For SARS-CoV-2, metadata columns must match the example exactly due to using Nextstrain's data resources
    * For bacteria, 'region' and 'country' are expected, and you can specify an additional column if desired, e.g. 'host type' or 'division'
    * For bacteria, the ParSNP core genome alignment tool renames strains to the FASTA filename and appends '.ref' to the reference genome. Metadata strain names must conform to this format.

### Example metadata files:

`metadata.txt` looks like (no header line):
| NC_045512v2 | REFERENCE |
| --- | --- |
| MT520428.1 | CONF_A |
| MT520429.1 | CONF_A |

`metadata_nextstrain.csv` looks like:
| strain | virus | date | region | country | division | 
| --- | --- | --- | --- | --- | --- |
| NC_045512v2 | ncov | 2019-12-22 | REFERENCE | China | ? |
| MT520428.1 | ncov | 2020-03-07 | North America | USA | CONF_A |
| MT520429.1 | ncov | 2020-03-08 | North America | USA | CONF_A |
