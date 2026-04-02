// codenarc-disable ModuleIncludedTwiceRule
include { WRITE_GFF3; WRITE_GFF3_BULK                                     } from "../../modules/output/gff3"
include { WRITE_JSON; WRITE_JSON_BULK                                     } from "../../modules/output/json"
include { WRITE_JSON as WRITE_JSONL; WRITE_JSON_BULK as WRITE_JSONL_BULK } from "../../modules/output/json"
include { WRITE_TSV; WRITE_TSV_BULK                                       } from "../../modules/output/tsv"
include { WRITE_XML; WRITE_XML_BULK                                       } from "../../modules/output/xml"

workflow OUTPUT {
    take:
    json            // [ meta, [json] ]
    seqdb
    formats
    outprefix
    nucleic
    interpro_version
    interproscan_version
    batch_size

    main:
    outfiles = Channel.empty()

    // convert to uppercase in case iprscan is imported directly into another workflow
    def formats_upper = formats.collect { it.toUpperCase() }

    if (batch_size < 50000) {
        if (formats_upper.contains("GFF3")) {
            WRITE_GFF3(json, file("${outprefix}.gff3"), seqdb, nucleic, interproscan_version)
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
    } else {
        if (formats_upper.contains("GFF3")) {
            WRITE_GFF3_BULK(json, file("${outprefix}.gff3"), seqdb, nucleic, interproscan_version)
            outfiles = outfiles.mix(WRITE_GFF3_BULK.out)
        }
        if (formats_upper.contains("JSON")) {
            WRITE_JSON_BULK(json, file("${outprefix}.json"), seqdb, nucleic, interproscan_version, interpro_version, false)
            outfiles = outfiles.mix(WRITE_JSON_BULK.out)
        }
        if (formats_upper.contains("JSONL")) {
            WRITE_JSONL_BULK(json, file("${outprefix}.jsonl"), seqdb, nucleic, interproscan_version, interpro_version, true)
            outfiles = outfiles.mix(WRITE_JSONL_BULK.out)
        }
        if (formats_upper.contains("TSV")) {
            WRITE_TSV_BULK(json, file("${outprefix}.tsv"), seqdb, nucleic)
            outfiles = outfiles.mix(WRITE_TSV_BULK.out)
        }
        if (formats_upper.contains("XML")) {
            WRITE_XML_BULK(json, file("${outprefix}.xml"), seqdb, nucleic, interproscan_version, interpro_version)
            outfiles = outfiles.mix(WRITE_XML_BULK.out)
        }
    }

    emit:
    outfiles.collect()
}
