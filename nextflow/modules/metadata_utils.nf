// Get metadata for clustertracker from Nextstrain metadata
process get_clustertracker_metadata_from_nextstrain {
    publishDir params.output_folder, mode: 'copy'
    container 'bash:devel-alpine3.17'

    cpus 1
    memory "1 GB"

    input:
        path nextstrain_metadata
        path clustertracker_focal_metadata

    output:
        path "clustertracker_metadata_with_context.tsv"

    shell:
        """
        #!/usr/bin/env bash
        set -eu
        # Extract strain and division columns from nextstrain metadata (cols 8, 9, 10 are region, country, division respectively) 
        # Remove header line
        awk -v OFS='\t' '{print \$1, \$10}' !{nextstrain_metadata} | tail -n +2 > clustertracker_metadata_from_nextstrain.tsv
        # Concatenate context metadata with focal metadata
        awk '{print}' !{clustertracker_focal_metadata} clustertracker_metadata_from_nextstrain.tsv > clustertracker_metadata_with_context.tsv
        """
}
