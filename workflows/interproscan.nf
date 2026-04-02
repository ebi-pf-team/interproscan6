// codenarc-disable ModuleIncludedTwiceRule
include { INIT_DATABASES       } from "../subworkflows/init/databases"
include { INIT_SEQUENCES       } from "../subworkflows/init/sequences"
include { LOOKUP               } from "../subworkflows/lookup"
include { SCAN_SEQUENCES as SCAN_REMAINING;
          SCAN_SEQUENCES       } from "../subworkflows/scan"
include { COMBINE              } from "../subworkflows/combine"
include { OUTPUT               } from "../subworkflows/output"

import uk.ac.ebi.interpro.InterProScan

workflow INTERPROSCAN {
    take:
    fasta_file              // Channel.fromPath(input fasta file)
    applications            // list[str], names of applications to run
    appl_config             // map, contents of the conf/applications.conf file
    data_dir                // path to the data directory
    outprefix               // path, prefix for output files
    formats                 // list[str], output formats
    interpro_version        // str, version of InterPro
    interproscan_version    // str, version of InterProScan
    interproscan_name       // str, name of the InterProScan workflow: "InterProScan6"
    no_matches_api          // boolean, use the Matches API
    matches_api_url         // str, url to the Matches API
    matches_api_chunk_size  // int, chunk size for Matches API requests
    matches_api_max_retries // int, max retries for Matches API requests
    batch_size              // int, number of sequences to process per batch
    sub_batch_size          // int, number of sequences per subbatch
    nucleic                 // boolean, input is nucleic acid sequences
    skip_repr_locations     // boolean, skip indentifying representative locations
    goterms                 // boolean, include GO terms in the output files
    pathways                // boolean, include pathway terms in the output files
    globus                  // boolean, use Globus to download databases
    enforce_compatibility   // boolean, whether to enforce compatibility between IPRScan and InterPro data versions

    main:

    INIT_DATABASES(
        applications,
        appl_config,
        data_dir,
        interpro_version,
        interproscan_version,
        !no_matches_api,
        globus,
        enforce_compatibility
    )
    appl_dirs        = INIT_DATABASES.out.directories
    interpro_version = INIT_DATABASES.out.interpro_version

    INIT_SEQUENCES(
        fasta_file,
        nucleic,
        batch_size
    )
    fasta = INIT_SEQUENCES.out.fasta
    seqdb = INIT_SEQUENCES.out.database

    scan_results = channel.empty()
    if (no_matches_api) {
        SCAN_SEQUENCES(
            fasta,
            appl_config,
            appl_dirs,
            applications,
            applications,
            sub_batch_size
        )
        scan_results = scan_results.mix(SCAN_SEQUENCES.out)
    }  else {
        LOOKUP(
            fasta,
            applications,
            matches_api_url,
            matches_api_chunk_size,
            matches_api_max_retries,
            interproscan_version
        )

        SCAN_REMAINING(
            LOOKUP.out.fasta,
            appl_config,
            appl_dirs,
            LOOKUP.out.appls_in_api.val,
            applications,
            sub_batch_size
        )

        SCAN_SEQUENCES(
            fasta,
            appl_config,
            appl_dirs,
            LOOKUP.out.appls_not_in_api.val,
            applications,
            sub_batch_size
        )

        scan_results = scan_results.mix(
            LOOKUP.out.json,
            SCAN_REMAINING.out,
            SCAN_SEQUENCES.out
        )
            .groupTuple()
            .map { meta, jsons -> meta
                tuple(meta, jsons.flatten())
            }
    }

    // Post-processing: aggregate results across all applications for each sequence and all cross-references
    COMBINE(
        scan_results,
        appl_config,
        appl_dirs,
        goterms,
        pathways,
        skip_repr_locations,
        batch_size
    )

    output_files = OUTPUT(
        COMBINE.out,
        seqdb,
        formats,
        outprefix,
        nucleic,
        interpro_version,
        interproscan_version,
        batch_size
    )

    emit:
    output_files
}

workflow PREPARE_INTERPROSCAN {
    // Parses an interproscan applications.config file and prepares it for use
    take:
    apps_config  // str repr of path to the applications.config file
    applications // list[str], names of applications to prepare
    use_gpu      // boolean, use GPU acceleration where applicable

    main:
    (apps_config, warn) = InterProScan.parseAppsConfig(use_gpu, applications, apps_config)
    if (warn) {
        log.warn warn
    }

    emit:
    apps_config
}
