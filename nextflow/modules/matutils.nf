nextflow.enable.dsl=2

// Align sequences to reference
process align_sequences {
    container 'pathogengenomics/usher@sha256:a311c7b896279a41a36608c29c8ecdc5b8420c01b27ffad160a1867098051c01'
    publishDir(path: "${params.output_folder}/matutils", mode: 'copy')

    cpus 2
    memory "1 GB"
    errorStrategy "retry"

    input:
        path input_fasta
        path reference

    output:
        path "alignment.fasta"

    shell:
        """
        set -eu
        mafft --thread 2 --auto --keeplength --addfragments !{input_fasta} !{reference} > alignment.fasta
        """
}

// Convert alignment to VCF file
process fasta_to_vcf {
    container 'pathogengenomics/usher@sha256:a311c7b896279a41a36608c29c8ecdc5b8420c01b27ffad160a1867098051c01'
    publishDir(path: "${params.output_folder}/matutils", mode: 'copy')

    cpus 1
    memory "1 GB"
    errorStrategy "retry"

    input:
        path alignment
        path reference
    
    output:
        path "masked.vcf"

    // TODO: abstract masking file
    shell:
        """
        set -eu
        wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf
        awk '{print}' !{reference} !{alignment} > alignment_w_reference.fasta
        faToVcf -maskSites=problematic_sites_sarsCov2.vcf alignment_w_reference.fasta masked.vcf
        """
}

// Build mutation-annotated tree, save in protobuf file
process build_mat {
    container 'pathogengenomics/usher@sha256:a311c7b896279a41a36608c29c8ecdc5b8420c01b27ffad160a1867098051c01'
    publishDir(path: "${params.output_folder}/matutils", mode: 'copy')

    cpus 2
    memory "1 GB"
    errorStrategy "retry"

    input:
        path vcf
        path tree

    output:
        path "mutation_annotated_tree.pb"

    shell:
        """
        set -eu
        usher -T 2 --vcf !{vcf} --tree !{tree} -o mutation_annotated_tree.pb
        """
}

// Apply heuristic to identify introductions
process matutils_introduce {
    container 'pathogengenomics/usher@sha256:a311c7b896279a41a36608c29c8ecdc5b8420c01b27ffad160a1867098051c01'
    publishDir(path: "${params.output_folder}/clustertracker", mode: 'copy')

    cpus 2
    memory "1 GB"
    errorStrategy "retry"

    input:
        path pb 
        path metadata 

    output:
        path "introductions.tsv"

    shell:
        """
        set -eu
        matUtils introduce --threads 2 -i !{pb} -s !{metadata} -o "introductions.tsv"
        """
}
