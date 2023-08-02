process mash_sketch {

    tag { output_basename }

    input:
      tuple val(sample_id), path(reads_1), path(reads_2)

    output:
      tuple val(sample_id), path ("${output_basename}_mash_sketch.txt")

    script:
      output_basename = reads_1.baseName.split('\\.')[0]
      """
      mash sketch -r ${reads_1} -m ${params.mash_sketch_minimum_copies} -k ${params.mash_sketch_kmer_size} -C ${output_basename} -o ${output_basename} &> ${output_basename}_mash_sketch.txt || true
      """
}

process mash_sketch_summary {

    tag { sample_id }

    cpus 1

    executor 'local'

    input:
      tuple val(sample_id), path(mash_sketch_output)

    output:
      tuple val(sample_id), path("${sample_id}_mash_sketch_summary.csv")

    script:
      """
      summarize_mash_sketch.py ${mash_sketch_output} --sample-id ${sample_id} > ${sample_id}_mash_sketch_summary.csv
      """
}
