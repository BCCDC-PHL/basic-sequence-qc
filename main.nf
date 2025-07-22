#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }                    from './modules/fastp.nf'
include { fastp_pre_dehosting }      from './modules/fastp.nf'
include { fastp_post_dehosting }     from './modules/fastp.nf'
include { dehost }                   from './modules/dehosting.nf'
include { combine_fastp_reports }    from './modules/dehosting.nf'


workflow {
  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq_input = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2']]] }
  } else {
    ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it.tail()] }
  }

  main:

    if (!params.dehost) {
        fastp(ch_fastq_input)
    }

    if (params.dehost) {
        fastp_pre_dehosting(ch_fastq_input)
	dehost(fastp_pre_dehosting.out.trimmed_reads)
	fastp_post_dehosting(dehost.out.dehosted_reads)
	combine_fastp_reports(fastp_pre_dehosting.out.metrics.join(fastp_post_dehosting.out.metrics))
    }

    output_prefix = params.prefix == '' ? params.prefix : params.prefix + '_'

    if (params.dehost) {
	combine_fastp_reports.out.metrics.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${output_prefix}basic_qc_stats.csv", storeDir: "${params.outdir}")
    } else {
	fastp.out.metrics.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${output_prefix}basic_qc_stats.csv", storeDir: "${params.outdir}")
    }

}
