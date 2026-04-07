include { WRITE_GFF3; WRITE_GFF3_BULK } from "${projectDir}/modules/output/gff3"
include { WRITE_JSON; WRITE_JSON_BULK } from "${projectDir}/modules/output/json"


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
    }
    
    emit:
    local
    bulk
}
