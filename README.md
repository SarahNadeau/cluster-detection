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
`metadata.txt` looks like (no header line):
| NC_045512v2 | REFERENCE |
| --- | --- |
| MT520428.1 | CONF_A |
| MT520429.1 | CONF_A |

`metadata_nextstrain.txt` looks like:
| strain | virus | date | region | country | division | 
| --- | --- | --- | --- | --- | --- |
| NC_045512v2 | ncov | 2019-12-22 | REFERENCE | China | ? |
| MT520428.1 | ncov | 2020-03-07 | North America | USA | CONF_A |
| MT520429.1 | ncov | 2020-03-08 | North America | USA | CONF_A |

Note that due to automated context set selection for SARS-CoV-2, the only fields available for trait reconstruction are 'region', 'country', and 'division'.
For the bacteria workflow, 'region' and 'country' are reconstructed by default and you can specify an additional column to reconstruct if you'd like, e.g. 'location' or 'annotation'.
The 'divison' metadata column is optional in this case.

Note also that the bacteria workflow uses parsnp, which renames strains to the filename and appends '.ref' to the reference filename. 
Your metadata strain names need to conform to the parsnp strain names.
