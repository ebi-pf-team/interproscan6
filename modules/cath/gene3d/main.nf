import groovy.json.JsonOutput

process SEARCH_GENE3D {
    label 'small', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), val(meta2), path(fasta)
    path hmmdir
    val hmmfile

    output:
    tuple val(meta), val(meta2), path("hmmsearch.out")

    script:
    """
    hmmsearch \
        -Z 65245 -E 0.001 \
        --cpu ${task.cpus} \
        ${hmmdir}/${hmmfile} ${fasta} > hmmsearch.out
    """
}

process RESOLVE_GENE3D {
    label 'tiny', 'ips6_container'

    input:
    tuple val(meta), val(meta2), path(hmmseach_out)

    output:
    tuple val(meta), val(meta2), path("resolved.out")

    script:
    """
    cath-resolve-hits \
        --input-format=hmmsearch_out \
        --min-dc-hmm-coverage=80 \
        --worst-permissible-bitscore=25 \
        --output-hmmer-aln ${hmmseach_out} > resolved.out
    """
}

process ASSIGN_CATH {
    label 'tiny', 'ips6_container'

    input:
    tuple val(meta), val(meta2), path(cath_resolve_out)
    path dirpath
    val dom2fam
    val disc_pickle

    output:
    tuple val(meta), val(meta2), path("cath.tsv")

    script:
    """
    python ${projectDir}/bin/cath/assign_cath_superfamilies.py \
        ${dirpath}/${dom2fam} \
        ${dirpath}/${disc_pickle} \
        ${cath_resolve_out} \
        cath.tsv
    """
}

process PARSE_CATHGENE3D {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(hmmseach_out), val(cath_tsv)

    output:
    tuple val(meta), val(meta2), path("cathgene3d.json")

    exec:
    def memberDb = "CATH-Gene3D"
    def hmmerMatches = HMMER3.parseOutput(hmmseach_out.toString(), memberDb)
    def cathDomains = CATH.parseAssignedFile(cath_tsv.toString())
    def matches = CATH.mergeWithHmmerMatches(cathDomains, hmmerMatches, memberDb)
    def outputFilePath = task.workDir.resolve("cathgene3d.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}