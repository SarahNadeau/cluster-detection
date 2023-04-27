nextflow.enable.dsl=2

// Align sequences to reference
process align_sequences {
    container 'pathogengenomics/usher@sha256:a311c7b896279a41a36608c29c8ecdc5b8420c01b27ffad160a1867098051c01'
    publishDir(path: "${params.output_folder}/matutils", mode: 'copy')

    label "process_low"

    input:
        path input_fasta
        path reference

    output:
        path "alignment.fasta"

    shell:
        """
        set -eu
        mafft --thread ${task.cpus} --auto --keeplength --addfragments !{input_fasta} !{reference} > alignment.fasta
        """
}

// Convert alignment to VCF file
process fasta_to_vcf {
    container 'pathogengenomics/usher@sha256:a311c7b896279a41a36608c29c8ecdc5b8420c01b27ffad160a1867098051c01'
    publishDir(path: "${params.output_folder}/matutils", mode: 'copy')

    input:
        path alignment
        path reference
        path mask_sites_vcf
    
    output:
        path "alignment.vcf"

    shell:
        """
        set -eu
        
        awk '{print}' !{reference} !{alignment} > alignment_w_reference.fasta

        # Convert to VCF, optionally masking some sites
        if [[ !{mask_sites_vcf} != 'NO_MASK_FILE' ]]; then
            faToVcf -maskSites=!{mask_sites_vcf} alignment_w_reference.fasta alignment.vcf
        else
            faToVcf alignment_w_reference.fasta alignment.vcf
        fi
        """
}

// Build mutation-annotated tree, save in protobuf file
process build_mat {
    container 'pathogengenomics/usher@sha256:a311c7b896279a41a36608c29c8ecdc5b8420c01b27ffad160a1867098051c01'
    publishDir(path: "${params.output_folder}/matutils", mode: 'copy')

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

    input:
        path pb 
        path metadata 

    output:
        path "introductions.tsv", emit: introductions_tsv
        path "clusters.tsv"

    shell:
        """
        set -eu
        matUtils introduce --threads 2 -i !{pb} -s !{metadata} -o "introductions.tsv" --cluster-output "clusters.tsv"
        """
}

// Convert mat protobuf file and metadata into taxonium format for interactive visualization
process pb_to_taxonium {
    container 'snads/taxoniumtools:2.0.91'
    publishDir(path: "${params.output_folder}/clustertracker", mode: 'copy')

    input:
        path pb
        path metadata

    output:
        path "taxonium.jsonl.gz"

    shell:
        """
        set -eu
        usher_to_taxonium \
            --input !{pb} \
            --output taxonium.jsonl.gz \
            --metadata !{metadata} \
            --key_column "sample" \
            --columns "introduction_node,introduction_rank,growth_score,cluster_size,intro_confidence,parent_confidence,region,origins,origins_confidence,mutation_path"
        """
}

// Extract introduction-descendent tips mapping from protobuf file
process pb_introductions_to_leaves {
    container 'snads/bte:0.9.0'
    publishDir(path: "${params.output_folder}/clustertracker", mode: 'copy')

    input:
        path pb
        path introductions

    output:
        path "node_to_leaves.csv"

    shell:
        template "get_node_to_leaves_mat.py"
}
