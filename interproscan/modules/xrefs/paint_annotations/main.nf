process PAINT_ANNOTATIONS {
    // Retrieve PAINT annotations for Panther hits
    // calculated and pre-calc becuase they are not retrieved from the Match Lookup
    label 'io'

    input:
    path matches
    path paint_anno_dir

    output:
    path "paint_anno_${matches}"

    script:
    """
    python3 $projectDir/interproscan/scripts/xrefs/paint_annotations.py ${matches} ${paint_anno_dir} > paint_anno_${matches}
    """
}
