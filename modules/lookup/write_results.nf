process WRITERESULTS {
    publishDir "${params.projectDir}/lookup_results", mode: 'copy'

    input:
    path 'parse_lookup'

    output:
    path '$projectDir/lookup'

    script:
    """
    if [ params.formats.contains('tsv')) ]; then
        python $projectDir/templates/write_output.py ${params.formats} -> $projectDir/lookup.tsv
    fi
    if [ params.formats.contains('json')) ]; then
        python $projectDir/templates/write_output.py ${params.formats} -> $projectDir/lookup.json
    fi
    if [ params.formats.contains('xml')) ]; then
        python $projectDir/templates/write_output.py ${params.formats} -> $projectDir/lookup.xml
    fi
    if [ params.formats.contains('gff3')) ]; then
        python $projectDir/templates/write_output.py ${params.formats} -> $projectDir/lookup.gff3
    fi
    """
}