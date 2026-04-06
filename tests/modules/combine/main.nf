include { COMBINE_MATCHES; COMBINE_MATCHES_BULK } from "${projectDir}/modules/combine"


workflow run_workflow {
    take: 
    json
    
    main:
    ch_default = COMBINE_MATCHES(json)
    ch_bulk    = COMBINE_MATCHES_BULK(json)
    
    ch_default.mix(ch_bulk).groupTuple(size: 2).view()

    emit:
    combined = ch_default.join(ch_bulk)
}
