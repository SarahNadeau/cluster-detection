// Build a phylogenetic tree using a maximum likelihood approach
process build_tree {
    container 'staphb/iqtree2:2.2.2.2'
    publishDir(path: "${params.output_folder}/iqtree", mode: 'copy')

    label "process_medium"

    input: 
        path alignment
        val outgroup_taxon

    output:
        path "tree.nwk"

    shell:
        """
        set -eu
        iqtree2 -s !{alignment} -o !{outgroup_taxon} -ntmax 2 -m HKY+F+G4 -pre 'tree'
        mv tree.treefile tree.nwk
        """
}