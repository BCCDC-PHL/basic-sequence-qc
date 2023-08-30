#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp } from './modules/fastp.nf'

workflow {
  ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it.tail()] }

  main:

    fastp(ch_fastq_input)

    output_prefix = params.prefix == '' ? params.prefix : params.prefix + '_'
    fastp.out.collectFile(keepHeader: true, sort: { it.text }, name: "${output_prefix}basic_qc_stats.csv", storeDir: "${params.outdir}")

}
