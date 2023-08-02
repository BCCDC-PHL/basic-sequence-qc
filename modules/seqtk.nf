process seqtk_fqchk {

    tag { sample_id }

    cpus 1

    input:
      tuple val(sample_id), path(reads_1), path(reads_2)

    output:
      tuple val(sample_id), path("${sample_id}_seqtk_fqchk_all_positions.csv"), path("${sample_id}_seqtk_fqchk_by_position.csv")

    script:
      """
      echo 'filename,position,num_bases,percent_a,percent_c,percent_g,percent_t,percent_n,average_q,error_q,percent_bases_below_q${params.seqtk_fqchk_threshold},percent_bases_above_q${params.seqtk_fqchk_threshold}' > header.csv
      seqtk fqchk -q ${params.seqtk_fqchk_threshold} ${reads_1} | tr \$'\\t' ',' | tail -n+3 | awk -F ',' 'BEGIN {OFS=FS}; {print "${reads_1}", \$0}' > ${sample_id}_R1_seqtk_fqchk_data.csv
      grep 'ALL' ${sample_id}_R1_seqtk_fqchk_data.csv > ${sample_id}_R1_seqtk_fqchk_data_all_positions.csv
      grep -v 'ALL' ${sample_id}_R1_seqtk_fqchk_data.csv > ${sample_id}_R1_seqtk_fqchk_data_by_position.csv || true
      seqtk fqchk -q ${params.seqtk_fqchk_threshold} ${reads_2} | tr \$'\\t' ',' | tail -n+3 | awk -F ',' 'BEGIN {OFS=FS}; {print "${reads_2}", \$0}' > ${sample_id}_R2_seqtk_fqchk_data.csv
      grep 'ALL' ${sample_id}_R2_seqtk_fqchk_data.csv > ${sample_id}_R2_seqtk_fqchk_data_all_positions.csv
      grep -v 'ALL' ${sample_id}_R2_seqtk_fqchk_data.csv > ${sample_id}_R2_seqtk_fqchk_data_by_position.csv || true
      cat header.csv ${sample_id}_R1_seqtk_fqchk_data_all_positions.csv ${sample_id}_R2_seqtk_fqchk_data_all_positions.csv > ${sample_id}_seqtk_fqchk_all_positions.csv
      cat header.csv ${sample_id}_R1_seqtk_fqchk_data_by_position.csv ${sample_id}_R2_seqtk_fqchk_data_by_position.csv > ${sample_id}_seqtk_fqchk_by_position.csv
      """
}

process seqtk_fqchk_summary {

    tag { sample_id }

    cpus 1

    executor 'local'

    input:
      tuple val(sample_id), path(seqtk_fqchk_output_all_positions), path(seqtk_fqchk_output_by_position)

    output:
      tuple val(sample_id), path("${sample_id}_seqtk_fqchk_summary.csv")

    script:
      """
      summarize_seqtk_fqchk.py ${seqtk_fqchk_output_all_positions} --sample-id ${sample_id} > ${sample_id}_seqtk_fqchk_summary.csv
      """
}
