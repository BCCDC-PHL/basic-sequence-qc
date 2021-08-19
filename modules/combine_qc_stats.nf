process combine_qc_stats {

  tag { sample_id }

  executor 'local'

  input:
    tuple val(sample_id), path(species_abundance), path(sequence_quality), path(estimated_coverage)

  output:
    tuple val(sample_id), path("${sample_id}_combined_qc_stats.csv")

  script:
  """
  printf "sample_id\\n${sample_id}\\n" > sample_id.csv
  awk -F ',' 'BEGIN {OFS=FS}; NR==1 { for (i=1; i<=NF; i++) {idx[\$i] = i} }; { print \$(idx["abundance_1_name"]), \$(idx["abundance_1_fraction_total_reads"]) }' ${species_abundance} | sed -s 's/abundance_1/most_abundant_species/g' > species_abundance_stats.csv
  awk -F ',' 'BEGIN {OFS=FS}; NR==1 { for (i=1; i<=NF; i++) {idx[\$i] = i} }; { print \$(idx["estimated_genome_size_bp"]), \$(idx["estimated_depth_coverage"]) }' ${estimated_coverage} > estimated_coverage_stats.csv
  awk -F ',' 'BEGIN {OFS=FS}; NR==1 { for (i=1; i<=NF; i++) {idx[\$i] = i} }; { print \$(idx["total_bases"]), \$(idx["average_base_quality"]), \$(idx["percent_bases_above_q${params.seqtk_fqchk_threshold}"]) }' ${sequence_quality} > sequence_quality_stats.csv
  paste -d ',' sample_id.csv species_abundance_stats.csv estimated_coverage_stats.csv sequence_quality_stats.csv > ${sample_id}_combined_qc_stats.csv
  """
}
