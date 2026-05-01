// codenarc-disable ModuleIncludedTwiceRule
include { WRITE_GFF3                } from "../../modules/write_output"
include { WRITE_JSON                } from "../../modules/write_output"
include { WRITE_JSON as WRITE_JSONL } from "../../modules/write_output"
include { WRITE_TSV                 } from "../../modules/write_output"
include { WRITE_XML                 } from "../../modules/write_output"

workflow OUTPUT {
    take:
    json            // [json]
    seqdb
    formats
    outprefix
    nucleic
    interpro_version
    interproscan_version

    main:
    outfiles = channel.empty()

    // convert to uppercase in case iprscan is imported directly into another workflow
    def formats_upper = formats.collect { it.toUpperCase() }

    if (formats_upper.contains("GFF3")) {
        WRITE_GFF3(json, file("${outprefix}.gff3"), seqdb, nucleic, interproscan_version, interpro_version)
        outfiles = outfiles.mix(WRITE_GFF3.out)
    }
    if (formats_upper.contains("JSON")) {
        WRITE_JSON(json, file("${outprefix}.json"), seqdb, nucleic, interproscan_version, interpro_version, false)
        outfiles = outfiles.mix(WRITE_JSON.out)
    }
    if (formats_upper.contains("JSONL")) {
        WRITE_JSONL(json, file("${outprefix}.jsonl"), seqdb, nucleic, interproscan_version, interpro_version, true)
        outfiles = outfiles.mix(WRITE_JSONL.out)
    }
    if (formats_upper.contains("TSV")) {
        WRITE_TSV(json, file("${outprefix}.tsv"), seqdb, nucleic)
        outfiles = outfiles.mix(WRITE_TSV.out)
    }
    if (formats_upper.contains("XML")) {
        WRITE_XML(json, file("${outprefix}.xml"), seqdb, nucleic, interproscan_version, interpro_version)
        outfiles = outfiles.mix(WRITE_XML.out)
    }

    emit:
    outfiles.collect()
}
