// Download Nextstrain-curated dataset
process download_nextstrain_covid_data {
    container 'nextstrain/base:build-20230411T103027Z'
    publishDir(
        path: "${params.output_folder}/augur", 
        mode: 'copy',
        pattern: 'metadata_*')

    cpus 1
    memory "1 GB"

    input:
        val region_name

    output:
        path "metadata_${region_name}.tsv.xz", emit: metadata
        path "aligned_${region_name}.fasta.xz", emit: alignment

    shell:
        """
        set -eu
        curl https://data.nextstrain.org/files/ncov/open/!{region_name}/metadata.tsv.xz -o metadata_!{region_name}.tsv.xz
        curl https://data.nextstrain.org/files/ncov/open/!{region_name}/aligned.fasta.xz -o aligned_!{region_name}.fasta.xz
        """
}

// Get proximity list of most related context sequences for a focal set
process get_proximities {
    container 'staphb/augur:16.0.3'
    publishDir(path: "${params.output_folder}/augur", mode: 'copy')

    cpus 2
    memory "1 GB"

    input: 
        path context_alignment
        path focal_alignment
        path reference
        val reference_name

    output:
        path "proximities.tsv"

    shell:
        """
        set -eu
        git clone https://github.com/nextstrain/ncov.git ./ncov
        /Python-3.8.0/python ./ncov/scripts/get_distance_to_focal_set.py \
            --alignment !{context_alignment} \
            --focal-alignment !{focal_alignment} \
            --reference !{reference} \
            --ignore-seqs !{reference_name} \
            --output proximities.tsv
        """
}

// Get priority ranking of most related context sequences for a focal set
process get_priorities {
    container 'staphb/augur:16.0.3'
    publishDir(path: "${params.output_folder}/augur", mode: 'copy')

    cpus 2
    memory "1 GB"

    input: 
        path context_alignment
        path proximities

    output:
        path "priorities.tsv", emit: priorities
        path "index.tsv", emit: index

    shell:
        """
        set -eu
        git clone https://github.com/nextstrain/ncov.git ./ncov

        augur index \
            --sequences !{context_alignment} \
            --output index.tsv

        /Python-3.8.0/python ./ncov/scripts/priorities.py \
            --sequence-index index.tsv \
            --proximities !{proximities} \
            --crowding-penalty 0 \
            --output priorities.tsv
        """
}

// Filter context sequence set based on region and genetic priorities
process augur_filter {
    container 'staphb/augur:16.0.3'
    publishDir(path: "${params.output_folder}/augur", mode: 'copy')

    cpus 2
    memory "1 GB"

    input: 
        path metadata
        path context_alignment
        path focal_alignment
        path priorities
        path index
        val min_sequence_length
        val group_by
        val max_context_sequences

    output:
        path "alignment_plus_filtered_context.fasta", emit: alignment_plus_filtered_context
        path "filtered_context_metadata.tsv", emit: filtered_context_metadata

    shell:
        """
        set -eu
        augur filter \
            --metadata !{metadata} \
            --sequences !{context_alignment} \
            --priority !{priorities} \
            --sequence-index !{index} \
            --min-length !{min_sequence_length} \
            --group-by !{group_by} \
            --subsample-max-sequences !{max_context_sequences} \
            --output filtered_context.fasta \
            --output-metadata filtered_context_metadata.tsv
        awk '{print}' !{focal_alignment} filtered_context.fasta > alignment_plus_filtered_context.fasta
        """
}
