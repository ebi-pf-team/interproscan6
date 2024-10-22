import groovy.json.JsonOutput

process RUN_SFLD {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out"), path("hmmsearch.dtbl"), path("hmmsearch.alignment")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -Z 378 --acc \
        --cut_ga \
        --cpu ${task.cpus} \
        -o hmmsearch.out \
        --domtblout hmmsearch.dtbl \
        -A hmmsearch.alignment \
        ${hmmdb} ${fasta}
    """
}

process POST_PROCESS_SFLD {
    label 'analysis_parser'

    input:
    tuple val(meta), val(hmmsearch_out), val(hmmsearch_dtbl), val(hmmsearch_alignment)
    val site_info

    output:
    tuple val(meta), path("sfld.post.processed.out")

    script:
    """
    $projectDir/interproscan/bin/sfld/sfld_postprocess \
        --alignment '${hmmsearch_alignment}' \
        --dom '${hmmsearch_dtbl}' \
        --hmmer-out '${hmmsearch_out}' \
        --site-info '${site_info}' \
        --output 'sfld.post.processed.out'
    """
}


process PARSE_SFLD {
    label 'analysis_parser'

    input:
    tuple val(meta), val(sfld_postprocess_output)
    val sfld_hierarchy_db

    output:
    tuple val(meta), val("sfld.json")

    exec:
    def outputFilePath = task.workDir.resolve("sfld.json")
    def matches = SFLD.parseOutput(sfld_postprocess_output.toString(), sfld_hierarchy_db.toString())
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
