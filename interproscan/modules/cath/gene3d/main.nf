import groovy.json.JsonOutput

process RESOLVE_GENE3D {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(hmmseach_out)

    output:
    tuple val(meta), path("resolved.out")

    script:
    """
    /opt/cath-tools/cath-resolve-hits \
        ${hmmseach_out} \
        --input-for hmmsearch_out \
        --min-dc-hmm-coverage=80 \
        --worst-permissible-bitscore 25 \
        --output-hmmer-aln > resolved.out
    """
}

process ASSIGN_CATH {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(cath_resolve_out)
    path dom2fam_file
    path disc_pickle_file

    output:
    tuple val(meta), path("cath.tsv")

    script:
    """
    python ${projectDir}/bin/cath/assign_cath_superfamilies.py \
        ${dom2fam_file} \
        ${disc_pickle_file} \
        ${cath_resolve_out} \
        cath.tsv
    """
}

process PARSE_CATHGENE3D {
    label 'run_locally'

    input:
    tuple val(meta), val(hmmseach_out), val(cath_tsv)

    output:
    tuple val(meta), path("cathgene3d.json")

    exec:
    def memberDb = "CATH-Gene3D"
    def hmmerMatches = HMMER3.parseOutput(hmmseach_out.toString(), memberDb)
    def cathDomains = CATH.parseAssignedFile(cath_tsv.toString())
    def matches = CATH.mergeWithHmmerMatches(cathDomains, hmmerMatches, memberDb)
    def outputFilePath = task.workDir.resolve("cathgene3d.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}