#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/deepvariant
========================================================================================
 nf-core/deepvariant Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/deepvariant
----------------------------------------------------------------------------------------
*/

nextflow.preview.dsl = 2

/*
 * Print help message if required
 */
if (params.help) {
    def command = "nextflow run nf-core/deepvariant --input samplesheet.csv -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$baseDir/nextflow_schema.json", command)
    exit 0
}

/*
 * Validate parameters
 */
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }

/*
 * Reference genomes
 */
params.fasta = params.genomes[params.genome]?.fasta
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) }

/*
 * Check parameters
 */
Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

/*
 * Print parameter summary
 */
// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    run_name = workflow.runName
}
summary = Schema.params_summary(workflow, params, run_name)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

workflow_summary = Schema.params_mqc_summary(summary)
ch_workflow_summary = Channel.value(workflow_summary)

/*
 * Include local pipeline modules
 */
include { OUTPUT_DOCUMENTATION } from './modules/local/output_documentation' params(params)
include { CHECK_SAMPLESHEET; check_samplesheet_paths } from './modules/local/check_samplesheet' params(params)


/*
 * Run the workflow
 */
workflow {

    CHECK_SAMPLESHEET(ch_input)
        .splitCsv(header:true, sep:',')
        .map { check_samplesheet_paths(it) }
        .set { ch_raw_reads }

}

/*
 * Send completion email
 */
workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary, run_name, baseDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}
