process add_reference {
    publishDir(path: "${params.output_folder}", mode: 'copy')

    input:
        path alignment
        path reference

    output:
        path "alignment.fasta"

    shell:
        """
        set -eu

        awk '{print}' !{reference} !{alignment} > alignment.fasta
        """
}

// Mask problematic sites in an alignment
// Using: https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/src/mask_alignment_using_vcf.py
// Credit to: https://virological.org/t/masking-strategies-for-sars-cov-2-alignments/480/2
process mask_alignment {
    publishDir(path: "${params.output_folder}", mode: 'copy')

    input:
        path alignment
        path mask_sites_vcf
        val reference_name

    output:
        path "masked_alignment.fasta"

    shell:
        """
        set -eu

        python ${projectDir}/../bin/mask_alignment_using_vcf.py \
            -i !{alignment} \
            -o masked_alignment.fasta \
            -v !{mask_sites_vcf} \
            --reference_id !{reference_name}
        """
}