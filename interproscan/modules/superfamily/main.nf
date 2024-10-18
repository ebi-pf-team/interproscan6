import groovy.json.JsonOutput

process SEARCH_SUPERFAMILY {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb
    path selfhits
    path cla
    path model
    path pdbj95d

    output:
    tuple val(meta), path("superfamily.out")

    script:
    """
    /opt/hmmer3/bin/hmmpress ${hmmdb}

    /opt/hmmer3/bin/hmmscan \
        -E 10 -Z 15438 \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmscan.out

    perl ${projectDir}/bin/superfamily/ass3_single_threaded.pl \
        -e 0.0001 -t n -f 1 \
        -s ${selfhits} \
        -r ${cla} \
        -m ${model} \
        -p ${pdbj95d} \
        ${fasta} \
        hmmscan.out \
        superfamily.out
    """
}