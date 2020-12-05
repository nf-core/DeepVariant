nextflow.enable.dsl = 2


// References
REF_FTPDIR = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids"
ref_fasta_gz = "${REF_FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ref_fasta_fai = "${REF_FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"


// Genomes
READS_FTPDIR = "ftp://ftp.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020"
BED_FILE = "${READS_FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.bed"
VCF_FILE = "${READS_FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz"
VCF_INDEX_FILE = "${READS_FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi"

// HG002 chr20 BAM
HTTPDIR = "https://storage.googleapis.com/deepvariant/case-study-testdata"
CHR_BAM = "${HTTPDIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
CHR_BAM_BAI = "${HTTPDIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai"

process DEEP_VARIANT {
    container "google/deepvariant:1.0.0"

    shell:

    '''

mkdir -p reference

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > reference/GRCh38_no_alt_analysis_set.fasta.fai


mkdir -p benchmark

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020

curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.bed > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.bed
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi


mkdir -p input
HTTPDIR=https://storage.googleapis.com/deepvariant/case-study-testdata

curl ${HTTPDIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam > input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam
curl ${HTTPDIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai > input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai


    run_deepvariant \
      --model_type WGS \
      --ref /reference/GRCh38_no_alt_analysis_set.fasta \
      --reads /input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
      --output_vcf /output/HG002.output.vcf.gz \
      --output_gvcf /output/HG002.output.g.vcf.gz \
      --num_shards $(nproc) \
      --regions chr20 
      
   '''

}


workflow {
    DEEP_VARIANT()
}
