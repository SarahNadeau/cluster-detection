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
        path "annotated_tree.nexus"
        path "confidence.csv"

    shell:
        """
        set -eu
        treetime mugration --tree !{tree} --states !{metadata} --attribute !{trait_name} --confidence --outdir .
        """
}
