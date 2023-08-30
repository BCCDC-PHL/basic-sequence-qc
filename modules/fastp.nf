process fastp {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}_fastp.csv")
    

    script:
    """
    fastp \
      -t ${task.cpus} \
      -i ${reads[0]} \
      -I ${reads[1]} \
      --cut_tail \
      -j ${sample_id}_fastp.json

    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
    """
}