// This workflow adds context sequences to a focal set of SARS-CoV-2 sequences using Nextstrain's curated data and augur filter

// Example parameters
params.input_fasta = "assets/test_data/sars_cov_2_boston_confa_outbreak_alignment.fasta"  // focal sequences, assumed aligned to same reference as Nextstrain context
params.context_region_name = "north-america"  // region to choose curated subsample (~4000 sequences) from; one of 'global', 'africa', 'asia', 'europe', 'north-america', 'oceania' or 'south-america'
params.reference_fasta = "assets/NC_045512v2.fa"  // reference that context and focal sequences are aligned to
params.reference_name  = "NC_045512v2"  // name of reference sequence in focal sequences to ignore for priority calculation
params.group_by = "'year', 'month', 'division'"  // string giving nextstrain augur filter group by specification (in the US divisions are states), e.g. "'year', 'month', 'division'"
params.sequences_per_group = 10  // number of sequences per group to filter to
params.output_folder = "results"  // where results files will be saved to

// Import processes from modules
include { download_nextstrain_covid_data; get_proximities; get_priorities; augur_filter } from './modules/augur.nf'

// The workflow itself
workflow {

    input_fasta = file(params.input_fasta)
    reference_fasta = file(params.reference_fasta)

    download_nextstrain_covid_data(params.context_region_name) 
    proximities = get_proximities(download_nextstrain_covid_data.out.alignment, input_fasta, reference_fasta, params.reference_name)
    get_priorities(download_nextstrain_covid_data.out.alignment, proximities)
    augur_filter(
        download_nextstrain_covid_data.out.metadata,
        download_nextstrain_covid_data.out.alignment,
        input_fasta,
        get_priorities.out.priorities, 
        get_priorities.out.index,
        params.group_by,
        params.sequences_per_group)

}
