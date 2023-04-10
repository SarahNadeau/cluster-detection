Data ideas:
* Aromatherapy outbreak where initially Sri Lanka was suspected as source, but was really India (sampling bias): https://www.ncbi.nlm.nih.gov/bioproject/PRJNA763213
* ClusterPicker data sets? No ground truth but all publicly available and include HIV, HCV, influenza
* NZ SARS-CoV-2 data annotated with community vs. border quarantine facility

Notes on files:

12864_2022_8936_MOESM4_ESM.ods is supplemental data from https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08936-4#additional-information where strain IDs (asked corresponding author for GISAID ID mapping) are annotated with known epidemiological cluster IDs.

sars-cov-2-SNF-A.tsv is metadata for a benchmark dataset originally from https://www.science.org/doi/10.1126/science.abe3261 and provided by Lee Katz. Based on the annotation "SNF A" and the code for the paper here https://github.com/JacobLemieux/sarscov2pub/blob/main/scripts/main_figures.R#L146 I assume these 63 samples are from the skilled nursing facility outbreak described in the paper. Therefore the 3 outbreak annotations are inferred from molecular clusters and are not the epidemiological ground truth.

nextstrain_sars-cov-2_boston_genbank_metadata.csv is downloaded from the Nextstrain build page https://auspice.broadinstitute.org/sars-cov-2/boston/genbank. It includes a column "CLF_A_Exposure"  with 63 entries which I assume corresponds to the same skilled nursing facility as above and a column "CONF_A_Exposure" with 28 entries which I assume correspond to the business conference outbreak described in the associated publication https://www.science.org/doi/10.1126/science.abe3261.

sars_cov_2_boston_snfa_outbreak.fasta I downloaded by pasting NCBI accession.version numbers from sars-cov-2-SNF-A.tsv into the nucleotide database search bar and sending results to a FASTA file.

sars_cov_2_boston_snfa_outbreak.fasta I downloaded the same as for the snfa outbreak except I first manually looked up each accession.version number on genbank based on the strain ID given in the Nextstrain metadata nextstrain_sars-cov-2_boston_genbank_metadata.csv.