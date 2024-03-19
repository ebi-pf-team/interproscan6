nextflow.enable.dsl=2
import org.yaml.snakeyaml.Yaml

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSE_SEQUENCE } from "$projectDir/modules/local/parse_sequence/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def f = file('input.yaml', checkifexists: true)
def input_yaml = new Yaml().load(f)

workflow {
    Channel.fromPath( input_yaml.input )
    .unique()
    .splitFasta( by: params.batchsize, file: true )
    .set { fasta_channel }

    PARSE_SEQUENCE(fasta_channel)
}
