// Align bacterial core genomes
process get_snps_and_tree {
    container 'staphb/parsnp:1.5.6'
    publishDir(path: "${params.output_folder}/parsnp", mode: 'copy')

    cpus 4
    memory "5 GB"

    input:
        path input_fasta_dir
        path reference

    output:
        path "parsnp/parsnp.vcf", emit: vcf
        path "parsnp/parsnp.tree", emit: tree
        path "snp_alignment.fasta", emit: snp_alignment

    shell:
        """
        set -eu

        # Align genomes, build tree (uses RAxML by default)
        parsnp -r !{reference} --sequences !{input_fasta_dir} --threads !{task.cpus} --output-dir parsnp --vcf

        # Extract SNP alignment from harvest file
        harvesttools -i parsnp/parsnp.ggr -S snp_alignment.fasta
        """

    stub:
        """
        mkdir parsnp
        cp ~/cluster-detection/nextflow/cached_results/klebsiella_urmc/parsnp/* parsnp/
        cp ~/cluster-detection/nextflow/cached_results/klebsiella_urmc/snp_alignment.fasta .
        """
}
