include { WRITE_GFF3; WRITE_GFF3_LOCAL                                     } from "../../modules/output/gff3"
include { WRITE_JSON; WRITE_JSON_LOCAL                                     } from "../../modules/output/json"
include { WRITE_JSON as WRITE_JSONL; WRITE_JSON_LOCAL as WRITE_JSONL_LOCAL } from "../../modules/output/json"
include { WRITE_TSV; WRITE_TSV_LOCAL                                       } from "../../modules/output/tsv"
include { WRITE_XML; WRITE_XML_LOCAL                                       } from "../../modules/output/xml"

workflow OUTPUT {
    take:
    ch_results
    seq_db_path
    formats
    outprefix
    nucleic
    iprscan_version
    db_releases
    batch_size

    main:
    outfiles = Channel.empty()

    // convert to uppercase in case iprscan is imported directly into another workflow
    def formats_upper = formats.collect { it.toUpperCase() }

    if (batch_size < 50000) {
        if (formats_upper.contains("GFF3")) {
            WRITE_GFF3_LOCAL(ch_results, "${outprefix}.gff3", seq_db_path, nucleic, iprscan_version)
            outfiles = outfiles.mix(WRITE_GFF3_LOCAL.out)
        }
        if (formats_upper.contains("JSON")) {
            WRITE_JSON_LOCAL(ch_results, "${outprefix}.json", seq_db_path, nucleic, iprscan_version, db_releases, false)
            outfiles = outfiles.mix(WRITE_JSON_LOCAL.out)
        }
        if (formats_upper.contains("JSONL")) {
            WRITE_JSONL_LOCAL(ch_results, "${outprefix}.jsonl", seq_db_path, nucleic, iprscan_version, db_releases, true)
            outfiles = outfiles.mix(WRITE_JSONL_LOCAL.out)
        }
        if (formats_upper.contains("TSV")) {
            WRITE_TSV_LOCAL(ch_results, "${outprefix}.tsv", seq_db_path, nucleic)
            outfiles = outfiles.mix(WRITE_TSV_LOCAL.out)
        }
        if (formats_upper.contains("XML")) {
            WRITE_XML_LOCAL(ch_results, "${outprefix}.xml", seq_db_path, nucleic, iprscan_version, db_releases)
            outfiles = outfiles.mix(WRITE_XML_LOCAL.out)
        }
    } else {
        if (formats_upper.contains("GFF3")) {
            WRITE_GFF3(ch_results, "${outprefix}.gff3", seq_db_path, nucleic, iprscan_version)
            outfiles = outfiles.mix(WRITE_GFF3.out)
        }
        if (formats_upper.contains("JSON")) {
            WRITE_JSON(ch_results, "${outprefix}.json", seq_db_path, nucleic, iprscan_version, db_releases, false)
            outfiles = outfiles.mix(WRITE_JSON.out)
        }
        if (formats_upper.contains("JSONL")) {
            WRITE_JSONL(ch_results, "${outprefix}.jsonl", seq_db_path, nucleic, iprscan_version, db_releases, true)
            outfiles = outfiles.mix(WRITE_JSONL.out)
        }
        if (formats_upper.contains("TSV")) {
            WRITE_TSV(ch_results, "${outprefix}.tsv", seq_db_path, nucleic)
            outfiles = outfiles.mix(WRITE_TSV.out)
        }
        if (formats_upper.contains("XML")) {
            WRITE_XML(ch_results, "${outprefix}.xml", seq_db_path, nucleic, iprscan_version, db_releases)
            outfiles = outfiles.mix(WRITE_XML.out)
        }
    }

    emit:
    outfiles.collect()
}
