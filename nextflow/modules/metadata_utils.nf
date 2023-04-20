// Get metadata (cols for strain, trait) from Nextstrain metadata
process get_metadata_from_nextstrain {
    publishDir(path: "${params.output_folder}", mode: 'copy')
    container 'nextstrain/base:build-20230411T103027Z'  // use any old container with /bin/bash

    cpus 1
    memory "1 GB"

    input:
        path nextstrain_metadata
        path focal_metadata

    output:
        path "metadata_with_context.tsv"

    shell:
        """
        set -eu

        # Extract strain and division columns from nextstrain metadata (cols 8, 9, 10 are region, country, division respectively)
        awk -v OFS='\t' 'BEGIN {FS = "\t"} ; {print \$1, \$10}' !{nextstrain_metadata} > metadata_from_nextstrain.tsv
        # Concatenate context metadata with focal metadata
        awk '{print}' metadata_from_nextstrain.tsv !{focal_metadata} > metadata_with_context_tmp.tsv
        # Remove any parentheses that will mess up tree file parsing
        sed 's/[()]//g' metadata_with_context_tmp.tsv > metadata_with_context.tsv
        """
}

// Concatenate nextstrain-formatted focal metadata with nextstrain context metadata
process add_metadata_to_nextstrain {
    publishDir(path: "${params.output_folder}", mode: 'copy')
    container 'nextstrain/base:build-20230411T103027Z'

    cpus 1
    memory "1 GB"

    input:
        path nextstrain_metadata
        path focal_metadata

    output:
        path "nextstrain_metadata_with_context.tsv"

    shell:
        """
        set -eu

        # Extract strain, virus, date, region, country, division columns from nextstrain metadata (cols 8, 9 are region, country respectively)
        awk -v OFS='\t' 'BEGIN {FS = "\t"} ; {print \$1, \$2, \$7, \$8, \$9, \$10}' !{nextstrain_metadata} | tail -n +2 > metadata_from_nextstrain.tsv
        # Concatenate context metadata with focal metadata
        awk '{print}' !{focal_metadata} metadata_from_nextstrain.tsv > nextstrain_metadata_with_context_tmp.tsv
        # Remove any parentheses that will mess up tree file parsing
        sed 's/[()]//g' nextstrain_metadata_with_context_tmp.tsv > nextstrain_metadata_with_context.tsv
        """
}

// Copy metadata to results folder
process save_metadata {
    publishDir(path: "${params.output_folder}", mode: 'copy')
    
    input:
        path metadata

    output:
        path ${metadata}

    shell:
        """
        echo "Copying metadata to results directory"
        """
}