include { COMBINE_MATCHES; COMBINE_MATCHES_BULK } from "${projectDir}/modules/combine"


workflow run_workflow {
    take: 
    json
    
    main:
    ch_default = COMBINE_MATCHES(json)
    ch_bulk    = COMBINE_MATCHES_BULK(json)
    
    emit:
    combined = ch_default.join(ch_bulk)
}
