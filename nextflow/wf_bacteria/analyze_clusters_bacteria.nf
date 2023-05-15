// This workflow is to analyze a bacterial outbreak dataset with several clustering methods.

// Default parameters (to run small test data)
params.input_fasta_dir = "test_data/samples"  // directory of *only* unaligned focal sequences
params.input_metadata = "test_data/metadata.txt"  // focal sequence metadata in clustertracker format (strain & location)
params.input_metadata_nextstrain = "test_data/metadata_nextstrain.csv"
params.reference_fasta = "test_data/reference.fasta"
params.trait_name = "location" // colname in metadata for trait to reconstruct
params.output_folder = "results_$workflow.start"  // where results files will be saved to
params.tn93_distance_threshold = 0.1 // genetic distance (under TN93 model) cutoff for clustering sequences (units are substitutions/site) 
params.hiv_trace_min_overlap = 1  // minimum number non-gap bases that must overlap for HIV-TRACE to calculate genetic distance (must be non-zero)
params.nextstrain_refine_params = "--coalescent opt --root reference.fasta.ref" // how should refine create a rooted timetree? See https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/refine.html

// Import processes from modules
include { get_snps_and_tree } from '../modules/parsnp.nf'
include { save_metadata } from '../modules/metadata_utils.nf'
include { build_mat; matutils_introduce; pb_to_taxonium; pb_introductions_to_leaves } from '../modules/matutils.nf'
include { hiv_trace} from '../modules/hiv_trace.nf'
include { augur_refine; augur_traits; augur_export } from '../modules/augur.nf'

// The workflow itself
workflow {

    input_fasta_dir = channel.fromPath(params.input_fasta_dir)
    input_metadata = channel.fromPath(params.input_metadata)
    input_metadata_nextstrain = channel.fromPath(params.input_metadata_nextstrain)
    reference_fasta = channel.fromPath(params.reference_fasta)

    // Align focal sequences, get SNP alignment in fasta & VCF, phylogeny
    get_snps_and_tree(input_fasta_dir, reference_fasta)

    // Save metadata
    save_metadata(input_metadata_nextstrain)

    // Run clustertracker to estimate introductions
    mat_pb = build_mat(get_snps_and_tree.out.vcf, get_snps_and_tree.out.tree)
    matutils_introduce(mat_pb, input_metadata)
    pb_to_taxonium(
        mat_pb, 
        matutils_introduce.out.introductions_tsv)
    pb_introductions_to_leaves(
        mat_pb, 
        matutils_introduce.out.introductions_tsv)

    // Run a SNP-distance based clustering method
    hiv_trace(
        get_snps_and_tree.out.snp_alignment,
        reference_fasta,
        params.tn93_distance_threshold,
        params.hiv_trace_min_overlap)

    // Run nextstrain augur
    augur_refine(
        get_snps_and_tree.out.tree,
        input_metadata_nextstrain,
        get_snps_and_tree.out.snp_alignment,
	    reference_fasta,
        params.nextstrain_refine_params)
    augur_traits(
        augur_refine.out.tree,
        input_metadata_nextstrain,
        params.trait_name)
    augur_export(
        augur_refine.out.tree,
        input_metadata_nextstrain,
        augur_traits.out.traits,
        augur_refine.out.branch_lengths,
        params.trait_name)

}
