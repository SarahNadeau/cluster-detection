// Example parameters
params.input_fasta = "assets/test_data/sars_cov_2_boston_confa_outbreak.fasta"  // sequences to be analyzed (don't need to be aligned)
params.metadata = "assets/test_data/sars_cov_2_boston_confa_outbreak_clustertracker_metadata.txt"  // tsv with names of samples and associated regions in 1st and 2nd columns
params.reference = "assets/NC_045512v2.fa"  // reference genome to align samples to
params.output_folder = "results"  // where results files will be saved to
params.outgroup_taxon = "NC_045512v2" // single taxon or list of taxon names for IQ-TREE outgroup

// Import processes from modules
include { align_sequences; fasta_to_vcf; build_mat; matutils_introduce } from './modules/matutils.nf'
include { build_tree } from './modules/iqtree.nf'

// The workflow itself
workflow {

    input_fasta = file(params.input_fasta)
    reference = file(params.reference)
    metadata = file(params.metadata)

    alignment = align_sequences(input_fasta, reference) 
    vcf = fasta_to_vcf(alignment, reference)
    tree = build_tree(alignment, params.outgroup_taxon)
    mat_pb = build_mat(vcf, tree)
    matutils_introduce(mat_pb, metadata)

}
