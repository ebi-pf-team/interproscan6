process REPRESENTATIVE_LOCATIONS {
    label     'mem_medium', 'time_short'
    container 'interpro/groovy:4.0.27-1'

    input:
    tuple val(meta), path(json, arity: '1..*', name: 'inputs/*')

    output:
    tuple val(meta), path('outputs/*json', arity: '1..*')

    script:
    """
    select-repr-locations.groovy \
        /opt/interproscan6/lib \
        inputs \
        outputs

    chmod -R 777 outputs
    """
}

// process REPRESENTATIVE_LOCATIONS {
//     label    'mem_medium', 'time_short'
//     executor 'local'

//     input:
//     tuple val(meta), val(json)

//     output:
//     tuple val(meta), path('repr_locs/*json', arity: '1..*')

//     exec:
//     def outdir = task.workDir.resolve("repr_locs")
//     uk.ac.ebi.interpro.ProcessReprLocations.run(json, outdir)
// }
