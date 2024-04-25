process SFLD_POST_PROCESSER {
    input:
        path out_file
        path out_dtbl
        path alignment
        val postprocessing_params // contains [0] bin and [1] site_annotations file path
        val tsv_pro

    output:
        path "${out_file}.post_processed.out"
        path "${out_dtbl}.post_processed.dtbl"
        path alignment
        val postprocessing_params // need to make up input var nums in hmmer_parser module

    script:
        """
        ${postprocessing_params[0]} \
            --alignments '${alignment}' \
            --dom '${out_dtbl}' \
            --hmmer-out '${out_file}' \
            --site-info '${postprocessing_params[1]}' \
            --output '${tsv_pro ? "${out_file}.processed.out" : "${out_dtbl}.processed.dtbl"}'

        python3 $projectDir/scripts/members/sfld/sfld_process_post_processed.py \
            '${tsv_pro ? "${out_file}.processed.out" : "${out_dtbl}.processed.dtbl"}' \
            ${tsv_pro ? "-O '${out_file}'" : "-d '${out_dtbl}'"}

        touch ${out_file}.post_processed.out
        touch ${out_dtbl}.post_processed.dtbl
        """
}