#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/deepvariant
========================================================================================
    Github : https://github.com/nf-core/deepvariant
    Website: https://nf-co.re/deepvariant
    Slack  : https://nfcore.slack.com/channels/deepvariant
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { DEEPVARIANT } from './workflows/deepvariant'

//
// WORKFLOW: Run main nf-core/deepvariant analysis pipeline
//
workflow NFCORE_DEEPVARIANT {
    DEEPVARIANT ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_DEEPVARIANT ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
