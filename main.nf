#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { seqtk_fqchk } from './modules/seqtk.nf'
include { seqtk_fqchk_summary } from './modules/seqtk.nf'
include { kraken2 } from './modules/kraken2.nf'
include { bracken } from './modules/bracken.nf'
include { abundance_top_n } from './modules/bracken.nf'
include { mash_sketch } from './modules/mash.nf'
include { mash_sketch_summary } from './modules/mash.nf'
include { combine_qc_stats } from './modules/combine_qc_stats.nf'


workflow {
  ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).filter{ !( it[0] =~ /Undetermined/ ) }.map{ it -> [it[0].split('_')[0], it[1], it[2]] }
  ch_kraken2_db = Channel.fromPath(params.kraken2_db)
  ch_bracken_db = Channel.fromPath(params.bracken_db)

  main:

    seqtk_fqchk(ch_fastq_input)
    seqtk_fqchk_summary(seqtk_fqchk.out)

    mash_sketch(ch_fastq_input)
    mash_sketch_summary(mash_sketch.out)

    kraken2(ch_fastq_input.combine(ch_kraken2_db))
    bracken(kraken2.out.combine(ch_bracken_db))

    abundance_top_n(bracken.out)

    output_prefix = params.prefix == '' ? params.prefix : params.prefix + '_'
    combine_qc_stats(abundance_top_n.out.map{ it -> [it[0], it[1]] }.join(seqtk_fqchk_summary.out.map{ it -> [it[0], it[1]] }).join(mash_sketch_summary.out)).map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${output_prefix}basic_qc_stats.csv", storeDir: "${params.outdir}")

}
