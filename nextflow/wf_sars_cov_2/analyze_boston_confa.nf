// This workflow is to analyze the Bonston confA outbreak dataset with several methods.

// Input & output data parameters
params.input_fasta = "../../clean_data/sars_cov_2_lemiux_boston/outbreak_samples_confa.fasta"  // unaligned focal sequences
params.input_metadata = "../../clean_data/sars_cov_2_lemiux_boston/clustertracker_metadata.txt"  // focal sequence metadata in clustertracker format (strain & location)
params.input_metadata_nextstrain = "../../clean_data/sars_cov_2_lemiux_boston/nextstrain_metadata.txt"  // focal sequence metadata in nextstrain format (strain & location)
params.trait_name = "division" // colname in metadata for trait to reconstruct
params.output_folder = "../results/boston_confa"  // where results files will be saved to

// Context data parameters
params.context_region_name = "north-america"  // draw context set from ~4000 Nextstrain-curated north america sequences
params.reference_fasta = "../../clean_data/sars_cov_2_lemiux_boston/reference_NC_045512v2.fa"  // SARS-CoV-2 reference genome
params.reference_name  = "NC_045512v2"  // name of reference sequence to ignore for priority calculation
params.min_sequence_length = 20000  // minimum length of context sequences
params.context_group_by = "'year', 'month', 'division'"  // Nextstrain augur filter group by specification for context set (in the US divisions are states)
params.max_context_sequences_divisions = 1000 // maximum overall number of sequences (will sample probabalistically from groups)
params.max_context_sequences_countries = 1000 // maximum overall number of sequences (will sample probabalistically from groups)

// Tree-building parameters
params.outgroup_taxon = "NC_045512v2" // root tree using reference sequence as outgroup

// Method-specific parameters
params.mask_sites_vcf = "../../clean_data/sars_cov_2_lemiux_boston/problematic_sites_sarsCov2.vcf" // vcf file of problematic sites to mask
params.tn93_distance_threshold = 0.0000667 // genetic distance (under TN93 model) cutoff for clustering sequences (units are substitutions/site)
params.hiv_trace_min_overlap = 1  // minimum number non-gap bases that must overlap for HIV-TRACE to calculate genetic distance (must be non-zero)

// Import processes from modules
include { download_nextstrain_covid_data; get_proximities; get_priorities; run_nextstrain_all } from '../modules/augur.nf'
include { augur_filter as filter_divisions; augur_filter as filter_countries; augur_aggregate } from '../modules/augur.nf'
include { build_tree } from '../modules/iqtree.nf'
include { align_sequences; fasta_to_vcf; build_mat; matutils_introduce; pb_to_taxonium } from '../modules/matutils.nf'
include { get_metadata_from_nextstrain; add_metadata_to_nextstrain } from '../modules/metadata_utils.nf'
include { treetime_mugration } from '../modules/treetime.nf'
include { hiv_trace } from '../modules/hiv_trace.nf'

// The workflow itself
workflow {

    input_fasta = Channel.fromPath(params.input_fasta)
    input_metadata = Channel.fromPath(params.input_metadata)
    input_metadata_nextstrain = Channel.fromPath(params.input_metadata_nextstrain)
    reference_fasta = Channel.fromPath(params.reference_fasta)
    mask_sites_vcf = Channel.fromPath(params.mask_sites_vcf)

    // Align focal sequences
    focal_alignment = align_sequences(input_fasta, reference_fasta)

    // Select context sequences
    download_nextstrain_covid_data(params.context_region_name) 
    proximities = get_proximities(
        download_nextstrain_covid_data.out.alignment, 
        focal_alignment, 
        reference_fasta, 
        params.reference_name)
    get_priorities(
        download_nextstrain_covid_data.out.alignment, 
        proximities)
    filter_divisions(
        download_nextstrain_covid_data.out.metadata,
        download_nextstrain_covid_data.out.alignment,
        get_priorities.out.priorities, 
        get_priorities.out.index,
        params.min_sequence_length,
        Channel.value("region == 'North America'"),
        Channel.value('division'),
        params.max_context_sequences_divisions)
    filter_countries(
        download_nextstrain_covid_data.out.metadata,
        download_nextstrain_covid_data.out.alignment,
        get_priorities.out.priorities, 
        get_priorities.out.index,
        params.min_sequence_length,
        Channel.value("region != 'North America'"),
        Channel.value('country'),
        params.max_context_sequences_countries)
    augur_aggregate(
        focal_alignment,
        filter_divisions.out.filtered_context.concat(filter_countries.out.filtered_context),
        filter_divisions.out.filtered_context_metadata.concat(filter_countries.out.filtered_context_metadata))

    // Get clustertracker metadata
    full_metadata = get_metadata_from_nextstrain(
        augur_aggregate.out.filtered_context_metadata, 
        input_metadata)

    // Run clustertracker to estimate introductions
    vcf = fasta_to_vcf(
        augur_aggregate.out.alignment_plus_filtered_context, 
        reference_fasta,
        mask_sites_vcf)
    tree = build_tree(
        augur_aggregate.out.alignment_plus_filtered_context, 
        params.outgroup_taxon)
    mat_pb = build_mat(vcf, tree)
    matutils_introduce(mat_pb, full_metadata)
    pb_to_taxonium(
        mat_pb, 
        matutils_introduce.out.introductions_tsv)

    // Run nextstrain mugration to estimate ancestral locations
    treetime_mugration(
        tree,
        full_metadata,
        params.trait_name)

    // Run a SNP-distance based clustering method
    hiv_trace(
        augur_aggregate.out.alignment_plus_filtered_context,
        reference_fasta,
        params.tn93_distance_threshold,
        params.hiv_trace_min_overlap)

    // Run nextstrain mugration in the context of a nextstrain workflow
    full_metadata_nextstrain = add_metadata_to_nextstrain(
        augur_aggregate.out.filtered_context_metadata,
        input_metadata_nextstrain)
    run_nextstrain_all(
        full_metadata_nextstrain,
        augur_aggregate.out.alignment_plus_filtered_context,
        params.trait_name)

    // TODO: run BEAST1 DTA to estiamte ancestral locations

}
