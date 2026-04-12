process REPRESENTATIVE_LOCATIONS {
    label    'mem_medium', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(json)

    output:
    tuple val(meta), path('repr_locs/*json', arity: '1..*')

    exec:
    def outdir = task.workDir.resolve("repr_locs")
    uk.ac.ebi.interpro.ProcessReprLocations.run(json, outdir)
}

