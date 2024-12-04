import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PREPARE_FUNFAM {
    label 'small'

    input:
    tuple val(meta), val(cathgene3d_json)
    val root_dir

    output:
    tuple val(meta), val(funfams)
    
    exec:
    File jsonFile = new File(cathgene3d_json.toString())
    Path rootDir = file(root_dir)
    JsonSlurper jsonSlurper = new JsonSlurper()
    funfams = jsonSlurper.parse(jsonFile)
        .values()
        .collect{ jsonMatches ->
            jsonMatches
                .values()
                .collect { jsonMatch ->
                    Match match = Match.fromMap(jsonMatch)
                    String accession = match?.signature?.accession                       
                    assert accession != null
                    assert accession.startsWith("G3DSA:")
                    String cathId = accession.substring(6)
                    return cathId
                }
        }
        .flatten()
        .unique()
        .findAll { cathId ->
            String hmmPath = cathId.split("\\.").join(File.separator) + ".hmm"
            Path fullPath = rootDir.resolve(hmmPath)
            return file(fullPath.toString()).isFile()
        }
}

process SEARCH_FUNFAM {
    label 'medium', 'ips6_container'

    input:
    tuple val(meta), path(fasta), val(supfams)
    path root_dir

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    def commands = ""
    supfams.each { cathId -> 
        String hmmFilePath = cathId.split("\\.").join(File.separator) + ".hmm"
        String hmmPath = "${root_dir.toString()}/${hmmFilePath}"
        commands += "/opt/hmmer3/bin/hmmsearch"
        commands += " -Z 65245 --cut_tc"
        commands += " --cpu ${task.cpus}"
        commands += " ${hmmPath} ${fasta} >> hmmsearch.out\n"
    }

    """
    ${commands}
    """
}

process RESOLVE_FUNFAM {
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


process PARSE_FUNFAM {
    label 'small'

    input:
    tuple val(meta), val(hmmseach_out), val(resolved_tsv)

    output:
    tuple val(meta), path("cathfunfam.json")

    exec:
    def hmmerMatches = HMMER3.parseOutput(hmmseach_out.toString())
    def funfamDomains = CATH.parseResolvedFile(resolved_tsv.toString())
    def matches = CATH.mergeWithHmmerMatches(funfamDomains, hmmerMatches)
    def outputFilePath = task.workDir.resolve("cathfunfam.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}