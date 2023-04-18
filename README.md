# Molecular cluster detection project

This project aims to be a methods comparison and guide for molecular cluster detection for public health.

## Repository structure
The project includes several case-study datasets in [clean_data](./clean_data/) and workflows to run several cluster detection methods on these data in [nextflow](./nextflow/).

Re-usable nextflow processes are in [modules](./nextflow/modules/).

## Inputs
Each dataset needs sequence data from a focal outbreak, context sequences for comparison, an appropriate reference genome for alignment, and metadata files.

### A reference file:  
* reference_NC_045512v2.fa

### Focal outbreak sequences:
* sars_cov_2_boston_confa_outbreak.fasta

### Context sequences
These depend on the outbreak scenario and can be:
1) downloaded by the pipeline, e.g. in [analyze_boston_confa.nf](./nextflow/analyze_boston_confa.nf)
2) stored in a directory with the focal outbreak samples, e.g. `clean_data/all_samples_for_analysis/*.fna`

### Metadata
`*_clustertracker_metadata.txt` looks like (no header line):
| NC_045512v2 | REFERENCE |
| --- | --- |
| MT520428.1 | CONF_A |
| MT520429.1 | CONF_A |

`*_nextstrain_metadata.txt` looks like:
| strain | virus | date | region | country | division |
| --- | --- | --- | --- | --- | --- |
| NC_045512v2 | ncov | 12/22/19 | REFERENCE | China | ? |
| MT520428.1 | ncov | 3/7/20 | CONF_A | USA | Massachusetts |
| MT520429.1 | ncov | 3/8/20 | CONF_A | USA | Massachusetts |
