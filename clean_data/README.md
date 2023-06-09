# These are notes on wrangling publicly accessible data files

## Data ideas:
* Aromatherapy outbreak where initially Sri Lanka was suspected as source, but was really India (sampling bias): https://www.ncbi.nlm.nih.gov/bioproject/PRJNA763213
* ClusterPicker data sets? No ground truth but all publicly available and include HIV, HCV, influenza

## Notes on files:

12864_2022_8936_MOESM4_ESM.ods 
* is supplemental data from https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08936-4#additional-information where strain IDs (asked corresponding author for GISAID ID mapping) are annotated with known epidemiological cluster IDs.

### Lemiux et al SARS-CoV-2 Boston

sars-cov-2-SNF-A.tsv 
* is metadata for a benchmark dataset originally from https://www.science.org/doi/10.1126/science.abe3261 and provided by Lee Katz. Based on the annotation "SNF A" and the code for the paper here https://github.com/JacobLemieux/sarscov2pub/blob/main/scripts/main_figures.R#L146 I assume these 63 samples are from the skilled nursing facility outbreak described in the paper. Therefore the 3 outbreak annotations are inferred from molecular clusters and are not the epidemiological ground truth.

nextstrain_sars-cov-2_boston_genbank_metadata.csv 
* is downloaded from the Nextstrain build page https://auspice.broadinstitute.org/sars-cov-2/boston/genbank. It includes a column "CLF_A_Exposure"  with 63 entries which I assume corresponds to the same skilled nursing facility as above and a column "CONF_A_Exposure" with 28 entries which I assume correspond to the business conference outbreak described in the associated publication https://www.science.org/doi/10.1126/science.abe3261.

sars_cov_2_lemiux_boston/outbreak_samples_snfa.fasta
* I downloaded by pasting NCBI accession.version numbers from sars-cov-2-SNF-A.tsv into the nucleotide database search bar and sending results to a FASTA file.

sars_cov_2_lemiux_boston/outbreak_samples_confa.fasta 
* I downloaded the same as for the snfa outbreak except I first manually looked up each accession.version number on genbank based on the strain ID given in the Nextstrain metadata nextstrain_sars-cov-2_boston_genbank_metadata.csv.

### Malek et al Klebsiella aerogenes URMC

GCF_000215745.1_ASM21574v1_genomic.fna
* Genomic DNA assembly, downloaded from NCBI; the Klebsiella aerogenes reference genome according to the publication

PRJNA504784
* Genomic DNA assembly, downloaded from NCBI; the Klebsiella aerogenes outbreak samples sequenced in the publication. I moved the files to the clean_data directory.
* In clean_data, I took the Genbank prefix out of the fasta headers using `for FILE in *.fna; do sed 's/>.*URMC/>URMC/g' $FILE > $FILE.noprefix; done`

context_samples_110/
* 110 context sequences from Malek et al. paper downloaded from Genbank
* In clean_data, I took the Genbank prefix out using `for FILE in *.fna; do sed 's/>.*aerogenes />/g' $FILE > $FILE.noprefix; done`

aac.02577-18-sd006.xlsx
* metadata for the 110 context sequences used in Malek et al.
* for analysis I assume all samples are from the 21st century (20XX-XX-XX) because otherwise nextstrain refine has a really hard time.

klebsiella_file_to_sample_name_mapping.txt
* hand-made in excel based on the other metadata files
* relevant for ParSNP, which names taxa by filename