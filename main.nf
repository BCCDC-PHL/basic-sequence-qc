#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hash_files as hash_fastq_input }  from './modules/hash_files.nf'
include { hash_files as hash_fastq_output } from './modules/hash_files.nf'
include { fastp }                           from './modules/fastp.nf'
include { fastp_pre_dehosting }             from './modules/fastp.nf'
include { fastp_post_dehosting }            from './modules/fastp.nf'
include { dehost }                          from './modules/dehosting.nf'
include { combine_fastp_reports }           from './modules/dehosting.nf'
include { pipeline_provenance }             from './modules/provenance.nf'
include { collect_provenance }              from './modules/provenance.nf'


workflow {

  ch_workflow_metadata = Channel.value([
    workflow.sessionId,
    workflow.runName,
    workflow.manifest.name,
    workflow.manifest.version,
    workflow.start,
  ])

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq_input = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2']]] }
  } else {
    ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it.tail()] }
  }
  ch_hostile_cache_dir = Channel.fromPath(params.hostile_cache_dir)

  main:

    hash_fastq_input(ch_fastq_input.combine(Channel.of("fastq-input")))

    // Each branch emits the same three channel shapes: the metrics csv, the reads
    // carried forward, and its provenance as a list of one file per process.
    if (params.dehost) {
        fastp_pre_dehosting(ch_fastq_input)
        dehost(fastp_pre_dehosting.out.trimmed_reads.combine(ch_hostile_cache_dir))
        fastp_post_dehosting(dehost.out.dehosted_reads)
        combine_fastp_reports(fastp_pre_dehosting.out.metrics.join(fastp_post_dehosting.out.metrics))

        ch_metrics    = combine_fastp_reports.out.metrics
        ch_reads_out  = dehost.out.dehosted_reads
        ch_provenance_qc = fastp_pre_dehosting.out.provenance.join(dehost.out.provenance).join(fastp_post_dehosting.out.provenance).map{ it -> [it[0], [it[1], it[2], it[3]]] }
    } else {
        fastp(ch_fastq_input)

        ch_metrics    = fastp.out.metrics
        ch_reads_out  = fastp.out.trimmed_reads
        ch_provenance_qc = fastp.out.provenance.map{ it -> [it[0], [it[1]]] }
    }

    // The hash records what the pipeline produced, whether or not it is published.
    hash_fastq_output(ch_reads_out.combine(Channel.of("fastq-output")))

    output_prefix = params.prefix == '' ? params.prefix : params.prefix + '_'
    ch_metrics.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${output_prefix}basic_qc_stats.csv", storeDir: "${params.outdir}")

    // Pipeline Provenance

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    // Per-process Provenance
    // We build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]

    ch_provenance = ch_fastq_input.map{ it -> it[0] }
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_fastq_input.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_provenance_qc).map{ it -> [it[0], it[1] + it[2]] }
    ch_provenance = ch_provenance.join(hash_fastq_output.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
}
