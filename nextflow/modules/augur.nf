// Download Nextstrain-curated dataset
process download_nextstrain_covid_data {
    container 'nextstrain/base:build-20230411T103027Z'
    publishDir(
        path: "${params.output_folder}/augur", 
        mode: 'copy',
        pattern: 'metadata_*')

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

    label "process_medium"

    input: 
        path context_alignment
        path focal_alignment
        path reference

    output:
        path "proximities.tsv"

    script:
        """
        set -eu
        
        REF_NAME="\$(grep '^>' ${reference} | sed 's/>//')"
        
        git clone https://github.com/nextstrain/ncov.git ./ncov
        
        /Python-3.8.0/python ./ncov/scripts/get_distance_to_focal_set.py \
            --alignment ${context_alignment} \
            --focal-alignment ${focal_alignment} \
            --reference ${reference} \
            --ignore-seqs \${REF_NAME} \
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

    label "process_medium"

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
        cp ~/cluster-detection/nextflow/cached_results/boston_confa/index.tsv .
	"""
}

// Get list of headers to definitely exclude when filtering for context
process get_context_exclude_list {

    input:
        path input_fasta
        path reference
    
    output:
        path "exclude.txt"

    shell:
        """
        grep "^>" !{reference} | sed 's/>//' > exclude.txt
        grep "^>" !{input_fasta} | sed 's/>//' >> exclude.txt
        """
}

// Filter context sequence set based on region and genetic priorities
process augur_filter {
    container 'staphb/augur:16.0.3'
    publishDir(path: "${params.output_folder}/augur", mode: 'copy')

    label "process_medium"

    input: 
        path metadata
        path context_alignment
        path priorities
        path index
        path exclude_strains
        val use_priorities
        val other_params

    output:
        tuple val("seqs"), path("filtered_context.fasta"), emit: filtered_context
        tuple val("metadatas"), path("filtered_context_metadata.tsv"), emit: filtered_context_metadata
        path("filtered_context_strains.txt"), emit: filtered_strains

    shell:
        """
        set -eu
        # handle priorities or no priorities
        if [[ !{use_priorities} == 'true' ]]; then
            echo "using Nextstrain priorities"
            augur filter \
                --metadata !{metadata} \
                --sequences !{context_alignment} \
                --sequence-index !{index} \
                --priority !{priorities} \
                --exclude !{exclude_strains} \
                !{other_params} \
                --output filtered_context.fasta \
                --output-metadata filtered_context_metadata.tsv \
                --output-strains filtered_context_strains.txt
        else
            echo "not using Nextstrain priorities"
            augur filter \
                --metadata !{metadata} \
                --sequences !{context_alignment} \
                --sequence-index !{index} \
                --exclude !{exclude_strains} \
                !{other_params} \
                --output filtered_context.fasta \
                --output-metadata filtered_context_metadata.tsv \
                --output-strains filtered_context_strains.txt
        fi
        """
}

// Aggregate focal alignment and filtered context to a final alignment
process augur_aggregate_2_filters {

    input:
        path focal_alignment
        tuple val(file_type), path('filtered_context1.fasta'), path('filtered_context2.fasta')
        tuple val(file_type), path('filtered_context_metadata1.tsv'), path('filtered_context_metadata2.tsv')

    output:
        path "alignment_plus_filtered_context.fasta", emit: alignment_plus_filtered_context
        path "all_filtered_context_metadata.tsv", emit: filtered_context_metadata

    shell:
        """
        awk '{print}' !{focal_alignment} filtered_context*.fasta > alignment_plus_filtered_context.fasta
        awk '{print}' filtered_context_metadata*.tsv > all_filtered_context_metadata.tsv
        """
}

// run the whole nextstrain workflow *after filtering*, including export to auspice for visualization
process run_nextstrain_all_vcf {
    container 'nextstrain/base:build-20230411T103027Z'
    containerOptions = "--user root"
    publishDir(path: "${params.output_folder}/augur", mode: 'copy')

    label "proces_medium"

    input: 
        path metadata
        path alignment
	    path reference
        val trait_name
        val other_refine_params

    output:
        path "auspice.json", emit: auspice_json
        path "tree_raw.nwk"
        path "tree.nwk"

    shell:
        """
        set -eu

        augur tree \
            --alignment !{alignment} \
	        --vcf-reference !{reference} \
            --output tree_raw.nwk

        augur refine \
            --tree tree_raw.nwk \
            --alignment !{alignment} \
	        --vcf-reference !{reference} \
            --metadata !{metadata} \
            --output-tree tree.nwk \
            --output-node-data branch_lengths.json \
            --timetree \
            !{other_refine_params}

        augur traits \
            --tree tree.nwk \
            --metadata !{metadata} \
            --output-node-data traits.json \
            --columns region country !{trait_name} \
            --confidence

        # skipping reconstruction of nucleotide and amino acid mutations

        augur export v2 \
            --tree tree.nwk \
            --metadata !{metadata} \
            --node-data branch_lengths.json \
                        traits.json \
            --color-by-metadata region country !{trait_name} \
            --output auspice.json
        """
}

process run_nextstrain_all {
    container 'nextstrain/base:build-20230411T103027Z'
    publishDir(path: "${params.output_folder}/augur", mode: 'copy')

    label "process_medium"

    input:
        path metadata
        path alignment
        val trait_name
        val other_refine_params

    output:
        path "auspice.json", emit: auspice_json
        path "tree_raw.nwk"
        path "tree.nwk"

    shell:
        """
        set -eu

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
            !{other_refine_params}

        augur traits \
            --tree tree.nwk \
            --metadata !{metadata} \
            --output-node-data traits.json \
            --columns region country !{trait_name} \
            --confidence

        # skipping reconstruction of nucleotide and amino acid mutations

        augur export v2 \
            --tree tree.nwk \
            --metadata !{metadata} \
            --node-data branch_lengths.json \
                        traits.json \
            --color-by-metadata region country !{trait_name} \
            --output auspice.json
        """
}
