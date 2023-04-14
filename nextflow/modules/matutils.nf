nextflow.enable.dsl=2

// Align sequences to reference
process align_sequences {
    container 'pathogengenomics/usher@sha256:a311c7b896279a41a36608c29c8ecdc5b8420c01b27ffad160a1867098051c01'
    publishDir(path: "${params.output_folder}/matutils", mode: 'copy')

    cpus 2
    memory "1 GB"

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

    input:
        path pb 
        path metadata 

    output:
        path "introductions.tsv", emit: introductions_tsv

    shell:
        """
        set -eu
        matUtils introduce --threads 2 -i !{pb} -s !{metadata} -o "introductions.tsv"
        """
}

// Convert mat protobuf file and metadata into taxonium format for interactive visualization
process pb_to_taxonium {
    container 'snads/taxoniumtools:2.0.91'
    publishDir(path: "${params.output_folder}/clustertracker", mode: 'copy')

    cpus 1
    memory "1 GB"

    input:
        path pb
        path metadata
        val key_colname
        val metadata_colnames

    output:
        path "taxonium.jsonl.gz"

    shell:
        """
        set -eu
        usher_to_taxonium \
            --input !{pb} \
            --output taxonium.jsonl.gz \
            --metadata !{metadata} \
            --key_column !{key_colname} \
            --columns !{metadata_colnames}
        """
}
