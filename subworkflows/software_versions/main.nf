workflow SOFTWARE_VERSIONS {
    take:
    versions
    pipeline_version
    commit_id
    nextflow_version
    outdir

    main:
    versions
        .mix(channel.of(
            tuple('Workflow', 'ebi-pf-team/interproscan6', uk.ac.ebi.interpro.SoftwareVersions.workflowVersion(pipeline_version, commit_id)),
            tuple('Workflow', 'Nextflow', nextflow_version.toString())
        ))
        .unique()
        .toList()
        .map { entries -> uk.ac.ebi.interpro.SoftwareVersions.toYaml(entries) }
        .collectFile(name: 'versions.yml', storeDir: outdir, sort: false)
        .set { report }

    emit:
    report
}
