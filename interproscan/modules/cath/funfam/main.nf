import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PREPARE_FUNFAM {
    input:
    tuple val(meta), val(cathgene3d_json)
    val root_dir 
    
    exec:
    File jsonFile = new File(cathgene3d_json.toString())
    Path rootDir = file(root_dir)
    JsonSlurper jsonSlurper = new JsonSlurper()
    def cathGene3dMatches = jsonSlurper.parse(jsonFile)
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
                    String hmmPath = cathId.split("\\.").join(File.separator) + ".hmm"
                    return rootDir.resolve(hmmPath)
                }
        }
        .flatten()
        .unique()
        .findAll { it.isFile() }

    println(cathGene3dMatches) 
}