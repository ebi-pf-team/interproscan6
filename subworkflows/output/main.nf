include { WRITE_JSON } from "../../modules/output/json"
include { WRITE_TSV  } from "../../modules/output/tsv"
include { WRITE_XML  } from "../../modules/output/xml"
include { WRITE_GFF3 } from "../../modules/output/gff3"

workflow OUTPUT {
    take:
    ch_results
    seq_db_path
    formats
    outprefix
    nucleic
    iprscan_version
    db_releases

    main:
    if (formats.contains("GFF3")) {
        WRITE_GFF3(ch_results, "${outprefix}.gff3", seq_db_path, nucleic, iprscan_version)
    }
    if (formats.contains("JSON")) {
        WRITE_JSON(ch_results, "${outprefix}.json", seq_db_path, nucleic, iprscan_version, db_releases)
    }
    if (formats.contains("TSV")) {
        WRITE_TSV(ch_results, "${outprefix}.tsv", seq_db_path, nucleic)
    }
    if (formats.contains("XML")) {
        WRITE_XML(ch_results, "${outprefix}.xml", seq_db_path, nucleic, iprscan_version, db_releases)
    }
}
