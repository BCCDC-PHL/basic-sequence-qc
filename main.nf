#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }                    from './modules/fastp.nf'
include { fastp as fastp_dehosted }  from './modules/fastp.nf'
include { dehost }                   from './modules/dehosting.nf'
include { combine_fastp_reports }    from './modules/dehosting.nf'


workflow {
  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq_input = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2']]] }
  } else {
    ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it.tail()] }
  }

  main:

    fastp(ch_fastq_input.combine(Channel.of("")))

    if (params.dehost) {
	dehost(ch_fastq_input)
	fastp_dehosted(dehost.out.dehosted_reads.combine(Channel.of("_dehosted")))
	combine_fastp_reports(fastp.out.metrics.join(fastp_dehosted.out.metrics))
    }

    output_prefix = params.prefix == '' ? params.prefix : params.prefix + '_'

    if (params.dehost) {
	combine_fastp_reports.out.metrics.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${output_prefix}basic_qc_stats.csv", storeDir: "${params.outdir}")
    } else {
	fastp.out.metrics.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${output_prefix}basic_qc_stats.csv", storeDir: "${params.outdir}")
    }

}
