// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options = initOptions(params.options)

process DEEP_VARIANT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
            mode: params.publish_dir_mode,
            saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: meta.id) }

    conda(params.enable_conda ? "bioconda::deepvariant=1.0.0" : null)
    container "google/deepvariant:1.0.0"

    input:
    tuple val(meta), val(sample_name), path(bam), path(bai)
    tuple val(regions), val(model_type)
    tuple path(ref_fasta), path(ref_fasta_fai)

    output:
    tuple val(meta), path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${sample_name}.g.vcf.gz"), path("${sample_name}.g.vcf.gz.tbi"), emit: gvcf
    tuple val(meta), path("*html"), emit: report

    script:

    """

    /opt/deepvariant/bin/run_deepvariant \
        --model_type ${model_type} \
        --ref ${ref_fasta} \
        --reads ${chr_bam} \
        --output_vcf ${sample_name}.vcf.gz \
        --output_gvcf ${sample_name}.g.vcf.gz \
        --num_shards ${task.cpus} \
        --regions ${regions}

        """
}
