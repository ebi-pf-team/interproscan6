include { WRITE_JSON } from "../../modules/output/json"
include { WRITE_TSV  } from "../../modules/output/tsv"
include { WRITE_XML  } from "../../modules/output/xml"

workflow OUTPUT {
    def fileName = params.input.split('/').last()
    def outFileName = "${params.outdir}/${fileName}"

    if (formats.contains("JSON")) {
        WRITE_JSON(ch_results, "${outFileName}", seq_db_path, params.nucleic, workflow.manifest.version)
    }
    if (formats.contains("TSV")) {
        WRITE_TSV(ch_results, "${outFileName}", seq_db_path, params.nucleic)
    }
    if (formats.contains("XML")) {
        WRITE_XML(ch_results, "${outFileName}", seq_db_path, params.nucleic, workflow.manifest.version)
    }
}