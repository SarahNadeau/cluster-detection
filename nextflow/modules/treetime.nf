// Infer ancestral trait values on a tree
// https://treetime.readthedocs.io/en/latest/tutorials/mugration.html
process treetime_mugration {
    container 'snads/treetime:0.9.4'
    publishDir(path: "${params.output_folder}/mugration", mode: 'copy')

    cpus 2
    memory "1 GB"

    input:
        path tree
        path metadata
        val trait_name
    
    output:
        path "annotated_tree.nexus", emit: annotated_tree_nexus
        path "confidence.csv"

    shell:
        """
        set -eu
        treetime mugration --tree !{tree} --states !{metadata} --attribute !{trait_name} --confidence --outdir .
        """
}

// Convert annotated tree output to NHX format for R
process convert_tree_to_nhx {
    publishDir(path: "${params.output_folder}/mugration", mode: 'copy')

    input:
        path annotated_tree_nexus

    output:
        path "annotated_tree.nhx"

    shell:
        """
        grep "^ Tree" !{annotated_tree_nexus} | sed "s/ Tree tree1=//g" | sed "s/&/&&NHX:/g" > annotated_tree_tmp.nhx
        # Edge case: get rid of parentheses that mess up parsing
        sed "s/Georgia (Asia)/Georgia Asia/g" annotated_tree_tmp.nhx > annotated_tree.nhx
        """
}
