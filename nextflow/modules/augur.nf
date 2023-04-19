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
        path "metadata_${region_name}.tsv.*z", emit: metadata
        path "aligned_${region_name}.fasta.xz", emit: alignment

    shell:
        """
        set -eu
        # handle special case where only .gz format metadata available for all GenBank
        if [[ !{region_name} == '.' ]]; then
            curl https://data.nextstrain.org/files/ncov/open/!{region_name}/metadata.tsv.gz -o metadata_!{region_name}.tsv.gz
        else
            curl https://data.nextstrain.org/files/ncov/open/!{region_name}/metadata.tsv.xz -o metadata_!{region_name}.tsv.xz
        fi
        curl https://data.nextstrain.org/files/ncov/open/!{region_name}/aligned.fasta.xz -o aligned_!{region_name}.fasta.xz
        """

    stub:
        """
        set -eu
        cp ~/cluster-detection/nextflow/cached_results/boston_confa/metadata_..tsv.gz .
        cp ~/cluster-detection/nextflow/cached_results/boston_confa/aligned_..fasta.xz .
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

    stub:
        """
        set -eu
        cp ~/cluster-detection/nextflow/cached_results/boston_confa/proximities.tsv .
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

    stub:
        """
        set -eu
        cp ~/cluster-detection/nextflow/cached_results/boston_confa/priorities.tsv .
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

// run the whole nextstrain workflow *after filtering*, including export to auspice for visualization
process run_nextstrain_all {
    container 'nextstrain/base:build-20230411T103027Z'
    publishDir(path: "${params.output_folder}/augur", mode: 'copy')

    cpus 2
    memory "1 GB"

    input: 
        path metadata
        path alignment

    output:
        path "auspice.json", emit: auspice_json
        path "tree_raw.nwk"
        path "tree.nwk"

    shell:
        """
        set -eu
        
        augur index \
            --sequences !{alignment} \
            --output sequence_index.tsv

        augur tree \
            --alignment !{alignment} \
            --output tree_raw.nwk

        augur refine \
            --tree tree_raw.nwk \
            --alignment !{alignment} \
            --metadata !{metadata} \
            --output-tree tree.nwk \
            --output-node-data branch_lengths.json \
            --timetree \
            --coalescent opt \
            --date-confidence \
            --date-inference marginal

        augur traits \
            --tree tree.nwk \
            --metadata !{metadata} \
            --output-node-data traits.json \
            --columns region country \
            --confidence

        # skipping reconstruction of nucleotide and amino acid mutations

        augur export v2 \
            --tree tree.nwk \
            --metadata !{metadata} \
            --node-data branch_lengths.json \
                        traits.json \
            --color-by-metadata region country \
            --output auspice.json
        """
}