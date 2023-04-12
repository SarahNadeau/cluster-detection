// Get metadata (cols for strain, trait) from Nextstrain metadata
process get_metadata_from_nextstrain {
    publishDir(path: "${params.output_folder}", mode: 'copy')
    container 'bash:devel-alpine3.17'

    cpus 1
    memory "1 GB"

    input:
        path nextstrain_metadata
        path focal_metadata

    output:
        path "metadata_with_context.tsv"

    shell:
        """
        #!/usr/bin/env bash
        set -eu

        # Extract strain and division columns from nextstrain metadata (cols 8, 9, 10 are region, country, division respectively)
        awk -v OFS='\t' '{print \$1, \$10}' !{nextstrain_metadata} > metadata_from_nextstrain.tsv
        # Concatenate context metadata with focal metadata
        awk '{print}' metadata_from_nextstrain.tsv !{focal_metadata} > metadata_with_context.tsv
        """
}
