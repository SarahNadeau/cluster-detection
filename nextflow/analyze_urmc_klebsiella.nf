// This workflow is to analyze the Bonston confA outbreak dataset with several methods.

// Input & output data parameters
params.input_fasta_dir = "../clean_data/klebsiella_malek_urmc/outbreak_samples_test"  // directory of *only* unaligned focal sequences
params.input_metadata = "../clean_data/klebsiella_malek_urmc/klebsiella_urmc_outbreak_metadata.txt"  // focal sequence metadata in clustertracker format (strain & location)
params.reference_fasta = "../clean_data/klebsiella_malek_urmc/reference_NC_015663v1.fna"
params.trait_name = "division" // colname in metadata for trait to reconstruct
params.output_folder = "results/urmc_klebsiella"  // where results files will be saved to

// Method-specific parameters
// NOTE: parsnp outputs SNP alignment, so the distance threshold (0.5% based on ad-hoc lit search) was corrected for approximate core-genome alignment length (4540000) and approximate SNP alignment length (86983) based on a test run
params.tn93_distance_threshold = 0.26 // genetic distance (under TN93 model) cutoff for clustering sequences (units are substitutions/site) 
params.hiv_trace_min_overlap = 1  // minimum number non-gap bases that must overlap for HIV-TRACE to calculate genetic distance (must be non-zero)

// Import processes from modules
include { get_snps_and_tree } from './modules/parsnp.nf'
include { build_mat; matutils_introduce } from './modules/matutils.nf'
include { treetime_mugration } from './modules/treetime.nf'
include { hiv_trace} from './modules/hiv_trace.nf'

// The workflow itself
workflow {

    input_fasta_dir = channel.fromPath(params.input_fasta_dir)
    input_metadata = channel.fromPath(params.input_metadata)
    reference_fasta = channel.fromPath(params.reference_fasta)

    // Align focal sequences, get SNP alignment in fasta & VCF, phylogeny
    get_snps_and_tree(input_fasta_dir, reference_fasta)

    // Run clustertracker to estimate introductions
    mat_pb = build_mat(get_snps_and_tree.out.vcf, get_snps_and_tree.out.tree)
    matutils_introduce(mat_pb, input_metadata)

    // Run nextstrain mugration to estimate ancestral locations
    treetime_mugration(
        get_snps_and_tree.out.tree,
        input_metadata,
        params.trait_name)

    // Run a SNP-distance based clustering method
    hiv_trace(
        get_snps_and_tree.out.snp_alignment,
        reference_fasta,
        params.tn93_distance_threshold,
        params.hiv_trace_min_overlap)

}