process COMBINE_MATCHES {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(json)

    output:
    tuple val(meta), path('combined/*.json', arity: '1..*')

    exec:
    def outdir = task.workDir.resolve("combined")
    uk.ac.ebi.interpro.ProcessCombine.run(json, outdir, 1000)
}
