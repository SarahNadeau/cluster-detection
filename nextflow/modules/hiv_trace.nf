// Cluster sequences by genetic distance under TN93 model
// Tool can generate an alignment when provided a reference, but here I assume sequences are already aligned
// Note: tool cannot handle arbitrary fasta headers, e.g. crashed with 'MT520287.1 Severe acute respiratory syndrome coronavirus 2 isolate SARS-CoV-2/human/USA/MA_MGH_00427/2020'
// Note: the 'fraction' argument does not apply, since I opt to average over all possible resolutions of ambiguities
// https://github.com/veg/hivtrace
process hiv_trace {
    container 'snads/hivtrace:0.5.0'
    publishDir(path: "${params.output_folder}/hiv_trace", mode: 'copy')

    label "proces_low"

    input:
        path input_fasta
        path reference
        val distance_threshold
        val min_overlap
    
    output:
        path "clusters.json"
        path "hivtrace.log"
        path "*.tn93output.csv"

    shell:
        """
        set -eu
        hivtrace \
            --input !{input_fasta} \
            --reference !{reference} \
            --threshold !{distance_threshold} \
            --ambiguities skip \
            --minoverlap !{min_overlap} \
            --fraction 0.25 \
            --skip-alignment \
            --output clusters.json
        """
}
