include { XREFS } from "$projectDir/interproscan/subworkflows/xrefs/xrefs.nf"

workflow XREFS {
    take:
    matches
    dataDir

    main:
    XREFS(matches, dataDir)

    emit:
    XREFS.out
}
