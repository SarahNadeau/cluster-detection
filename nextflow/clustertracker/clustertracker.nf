nextflow.enable.dsl=2

// Script parameters
params.pb = "/Users/snadeau/Documents/repos/cluster-detection/nextflow/clustertracker/assets/test_data/public-2021-06-09.all.masked.nextclade.pangolin.pb.gz"  // protobuf
params.metadata = "/Users/snadeau/Documents/repos/cluster-detection/nextflow/clustertracker/assets/test_data/regional-samples.txt"  // tsv with names of samples and associated regions in 1st and 2nd columns
params.output_folder = "results"

process matutils_introduce {
    publishDir params.output_folder

    container "pathogengenomics/usher:latest"
    cpus 1
    memory "1 GB"
    errorStrategy "retry"

    input:
        path pb 
        path metadata 

    output:
        path "test.tsv"

    shell:
        """
        set -eu

        echo `pwd`
        matUtils introduce -i !{pb} -s !{metadata} -o "test.tsv"
        """
}

workflow {
    pb_ch = file(params.pb)
    metadata_ch = file(params.metadata)

    matutils_introduce(pb_ch, metadata_ch)
}