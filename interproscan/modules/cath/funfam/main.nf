import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PREPARE_FUNFAM {
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
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta), val(supfams)
    path root_dir

    output:
    tuple val(meta), stdout

    script:
    def commands = ""
    supfams.each { cathId -> 
        String hmmFilePath = cathId.split("\\.").join(File.separator) + ".hmm"
        String hmmPath = "${root_dir.toString()}/${hmmFilePath}"
        commands += "/opt/hmmer3/bin/hmmsearch"
        commands += " -Z 65245 --cut_tc"
        commands += " --cpu ${task.cpus}"
        commands += " ${hmmPath} ${fasta}\n"
    }

    """
    ${commands}
    """
}