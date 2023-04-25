// This workflow is to analyze the URMC K. aerogenes outbreak dataset with several methods.

// Input & output data parameters
params.input_fasta_dir = "../../clean_data/klebsiella_malek_urmc/all_samples_for_analysis"  // directory of *only* unaligned focal sequences
params.input_metadata = "../../clean_data/klebsiella_malek_urmc/clustertracker_metadata.txt"  // focal sequence metadata in clustertracker format (strain & location)
params.input_metadata_nextstrain = "../../clean_data/klebsiella_malek_urmc/nextstrain_metadata.csv"
params.reference_fasta = "../../clean_data/klebsiella_malek_urmc/reference_NC_015663v1.fna"
params.trait_name = "location" // colname in metadata for trait to reconstruct
params.output_folder = "../results/urmc_klebsiella_$workflow.start"  // where results files will be saved to

// Method-specific parameters
// NOTE: parsnp outputs SNP alignment, so the distance threshold (0.5% based on ad-hoc lit search) was corrected for approximate core-genome alignment length (4540000) and approximate SNP alignment length (86983) based on a test run
params.tn93_distance_threshold = 0.1 // genetic distance (under TN93 model) cutoff for clustering sequences (units are substitutions/site) 
params.hiv_trace_min_overlap = 1  // minimum number non-gap bases that must overlap for HIV-TRACE to calculate genetic distance (must be non-zero)
params.reference_name  = "reference_NC_015663v1.fna.ref"  // name of reference sequence to root augur tree on (note parsnp appends ".ref")

// Import processes from modules
include { get_snps_and_tree } from '../modules/parsnp.nf'
include { save_metadata } from '../modules/metadata_utils.nf'
include { build_mat; matutils_introduce; pb_to_taxonium } from '../modules/matutils.nf'
include { treetime_mugration } from '../modules/treetime.nf'
include { hiv_trace} from '../modules/hiv_trace.nf'
include { run_nextstrain_all_vcf } from '../modules/augur.nf'

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

    // Run nextstrain mugration to estimate ancestral locations
    treetime_mugration(
        get_snps_and_tree.out.tree,
        input_metadata_nextstrain,
        params.trait_name)

    // Run a SNP-distance based clustering method
    hiv_trace(
        get_snps_and_tree.out.snp_alignment,
        reference_fasta,
        params.tn93_distance_threshold,
        params.hiv_trace_min_overlap)

    // Run nextstrain mugration in the context of a nextstrain workflow
    run_nextstrain_all_vcf(
        input_metadata_nextstrain,
        get_snps_and_tree.out.vcf,
	    reference_fasta,
        params.trait_name,
        Channel.value("--root ${params.reference_name}"))

}
