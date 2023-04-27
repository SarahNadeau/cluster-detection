// This workflow is to analyze a SARS-CoV-2 outbreak dataset with several clustering methods.

// Default parameters (to run small test data)
params.input_fasta = "test_data/samples.fasta"  // unaligned focal sequences
params.input_metadata = "test_data/metadata.txt"  // focal sequence metadata in clustertracker format (strain & location)
params.input_metadata_nextstrain = "test_data/metadata_nextstrain.csv"  // focal sequence metadata in nextstrain format (strain & location)
params.trait_name = "division" // colname in metadata for trait to reconstruct
params.output_folder = "test_results_$workflow.start"  // where results files will be saved to
params.context_region_name = "global"  // draw context set from ~4000 Nextstrain-curated north america sequences
params.reference_fasta = "test_data/reference.fasta"  // SARS-CoV-2 reference genome
params.reference_name  = "NC_045512v2"  // name of reference sequence to ignore for priority calculation
<<<<<<< HEAD:nextflow/wf_sars_cov_2/analyze_boston_confa.nf
params.exclude_strains = "../../clean_data/sars_cov_2_lemiux_boston/exclude_augur_filter.txt"  // exclude focal samples from context set selection
params.max_similarity_seqs = 200  // number most genetically similar sequences (from anywhere) to include in context
params.max_geocontext_seqs = 400  // number of geographic context sequences to include in context (will be divided by region, year, month)

// Tree-building parameters
=======
params.filter_similarity_specs = "--min-length 20000 --max-date '2020-07-01' --subsample-max-sequences 10"
params.filter_geocontext_specs = "--min-length 20000 --max-date '2020-07-01' --group-by region year month --subsample-max-sequences 20"
>>>>>>> traverse-mat:nextflow/wf_sars_cov_2/analyze_clusters_sars_cov_2.nf
params.outgroup_taxon = "NC_045512v2" // root tree using reference sequence as outgroup
params.mask_sites_vcf = "assets/problematic_sites_sarsCov2.vcf" // vcf file of problematic sites to mask
params.tn93_distance_threshold = 0.0000667 // genetic distance (under TN93 model) cutoff for clustering sequences (units are substitutions/site)
params.hiv_trace_min_overlap = 1  // minimum number non-gap bases that must overlap for HIV-TRACE to calculate genetic distance (must be non-zero)
params.nextstrain_refine_params = "--coalescent opt --date-confidence --keep-polytomies --clock-rate 0.0008 --clock-std-dev 0.004 --root NC_045512v2" // how should refine create a rooted timetree? See https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/refine.html

// Import processes from modules
include { download_nextstrain_covid_data; get_proximities; get_priorities; run_nextstrain_all } from '../modules/augur.nf'
include { augur_filter as filter_1; augur_filter as filter_2; augur_aggregate_2_filters; get_context_exclude_list } from '../modules/augur.nf'
include { build_tree } from '../modules/iqtree.nf'
include { align_sequences; fasta_to_vcf; build_mat; matutils_introduce; pb_to_taxonium; pb_introductions_to_leaves } from '../modules/matutils.nf'
include { get_metadata_from_nextstrain; add_metadata_to_nextstrain; combine_exclude_files } from '../modules/metadata_utils.nf'
include { treetime_mugration; convert_tree_to_nhx } from '../modules/treetime.nf'
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
        reference_fasta)
    get_priorities(
        download_nextstrain_covid_data.out.alignment, 
        proximities)
    exclude_strains = get_context_exclude_list(
        focal_alignment,
        reference_fasta)
    filter_1(
        download_nextstrain_covid_data.out.metadata,
        download_nextstrain_covid_data.out.alignment,
        get_priorities.out.priorities, 
        get_priorities.out.index,
        exclude_strains,
        Channel.value("true"), 
        params.filter_similarity_specs)
    combine_exclude_files(
        exclude_strains,
        filter_1.out.filtered_strains) // don't re-select similarity context in geographic context set
    filter_2(
        download_nextstrain_covid_data.out.metadata,
        download_nextstrain_covid_data.out.alignment,
        get_priorities.out.priorities, 
        get_priorities.out.index,
        combine_exclude_files.out.exclude_strains_combined,
        Channel.value("false"), // don't use genetic similarity here
        params.filter_geocontext_specs)

    augur_aggregate_2_filters(
        focal_alignment,
        filter_1.out.filtered_context.join(filter_2.out.filtered_context),
        filter_1.out.filtered_context_metadata.join(filter_2.out.filtered_context_metadata))

    // Get clustertracker metadata
    full_metadata = get_metadata_from_nextstrain(
        augur_aggregate_2_filters.out.filtered_context_metadata, 
        input_metadata)

    // Run clustertracker to estimate introductions
    vcf = fasta_to_vcf(
        augur_aggregate_2_filters.out.alignment_plus_filtered_context, 
        reference_fasta,
        mask_sites_vcf)
    tree = build_tree(
        augur_aggregate_2_filters.out.alignment_plus_filtered_context, 
        params.outgroup_taxon)
    mat_pb = build_mat(vcf, tree)
    matutils_introduce(mat_pb, full_metadata)
    pb_to_taxonium(
        mat_pb, 
        matutils_introduce.out.introductions_tsv)
    pb_introductions_to_leaves(
        mat_pb, 
        matutils_introduce.out.introductions_tsv)

    // Run nextstrain mugration to estimate ancestral locations
    treetime_mugration(
        tree,
        full_metadata,
        params.trait_name)
    convert_tree_to_nhx(treetime_mugration.out.annotated_tree_nexus)

    // Run a SNP-distance based clustering method
    hiv_trace(
        augur_aggregate_2_filters.out.alignment_plus_filtered_context,
        reference_fasta,
        params.tn93_distance_threshold,
        params.hiv_trace_min_overlap)

    // Run nextstrain mugration in the context of a nextstrain workflow
    full_metadata_nextstrain = add_metadata_to_nextstrain(
        augur_aggregate_2_filters.out.filtered_context_metadata,
        input_metadata_nextstrain)
    run_nextstrain_all(
        full_metadata_nextstrain,
        augur_aggregate_2_filters.out.alignment_plus_filtered_context,
        params.trait_name,
        params.nextstrain_refine_params)

    // TODO: run BEAST1 DTA to estiamte ancestral locations

}
