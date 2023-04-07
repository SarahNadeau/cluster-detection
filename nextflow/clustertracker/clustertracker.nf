nextflow.enable.dsl=2

// Script parameters
params.reference = "assets/NC_045512v2.fa"
params.metadata = "assets/test_data/regional-samples.txt"  // tsv with names of samples and associated regions in 1st and 2nd columns
params.output_folder = "results"

// Align sequences to reference
process align_sequences {
    publishDir params.output_folder

    cpus 2
    memory "1 GB"
    errorStrategy "retry"

    input:
        path input_fasta
        path reference

    output:
        path "alignment.fasta"

    shell:
        """
        set -eu
        mafft --thread 2 --auto --keeplength --addfragments !{input_fasta} !{reference} > alignment.fasta
        """
}

// process matutils_introduce {
//     publishDir params.output_folder

//     cpus 1
//     memory "1 GB"
//     errorStrategy "retry"

//     input:
//         path pb 
//         path metadata 

//     output:
//         path "test.tsv"

//     shell:
//         """
//         set -eu
//         matUtils introduce -i !{pb} -s !{metadata} -o "test.tsv"
//         """
// }

workflow {
    input_fasta_ch = file(params.input_fasta)
    reference_ch = file(params.reference)

    align_sequences(input_fasta_ch, reference_ch)

    // metadata_ch = file(params.metadata)
    // matutils_introduce(pb_ch, metadata_ch)
}