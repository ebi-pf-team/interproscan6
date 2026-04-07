include { WRITE_GFF3; WRITE_GFF3_BULK } from "${projectDir}/modules/output/gff3"
include { WRITE_JSON; WRITE_JSON_BULK } from "${projectDir}/modules/output/json"
include { WRITE_TSV; WRITE_TSV_BULK } from "${projectDir}/modules/output/tsv"
include { WRITE_XML; WRITE_XML_BULK   } from "${projectDir}/modules/output/xml"


workflow run_workflow {
    take: 
    format
    json
    output
    seqdb
    is_nucleic
    iprscan_version
    interpro_version
   
    main:
    if (format == "gff3") {
        local = WRITE_GFF3(json, output, seqdb, is_nucleic, iprscan_version)
        bulk  = WRITE_GFF3_BULK(json, output, seqdb, is_nucleic, iprscan_version)
    } else if (format == "json") {
        local = WRITE_JSON(json, output, seqdb, is_nucleic, iprscan_version, interpro_version, false)
        bulk  = WRITE_JSON_BULK(json, output, seqdb, is_nucleic, iprscan_version, interpro_version, false)
    } else if (format == "jsonl") {
        local = WRITE_JSON(json, output, seqdb, is_nucleic, iprscan_version, interpro_version, true)
        bulk  = WRITE_JSON_BULK(json, output, seqdb, is_nucleic, iprscan_version, interpro_version, true)        
    } else if (format == "tsv") {
        local = WRITE_TSV(json, output, seqdb, is_nucleic)
        bulk  = WRITE_TSV_BULK(json, output, seqdb, is_nucleic)        
    } else if (format == "xml") {
        local = WRITE_XML(json, output, seqdb, is_nucleic, iprscan_version, interpro_version)
        bulk  = WRITE_XML_BULK(json, output, seqdb, is_nucleic, iprscan_version, interpro_version)        
    }
    
    emit:
    local
    bulk
}
