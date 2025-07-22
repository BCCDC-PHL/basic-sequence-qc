process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*_fastp.*", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_trimmed_R*.fastq.gz", mode: 'copy', enabled: params.publish_trimmed_reads

    input:
    tuple val(sample_id), path(reads), val(suffix)

    output:
    tuple val(sample_id), path("${sample_id}${suffix}_fastp.csv")   , emit: metrics
    tuple val(sample_id), path("${sample_id}_fastp.json")           , emit: report_json
    tuple val(sample_id), path("${sample_id}_fastp.html")           , emit: report_html
    tuple val(sample_id), path("${sample_id}_trimmed_R*.fastq.gz")  , emit: trimmed_reads

    script:
    worker_threads = task.cpus - 1
    """
    fastp \
      --thread ${worker_threads} \
      -i ${reads[0]} \
      -I ${reads[1]} \
      --cut_tail \
      --trim_poly_g \
      --overrepresentation_analysis \
      --detect_adapter_for_pe \
      -o ${sample_id}_trimmed_R1.fastq.gz \
      -O ${sample_id}_trimmed_R2.fastq.gz \
      --report_title "fastp report: ${sample_id}" \
      --json ${sample_id}_fastp.json \
      --html ${sample_id}_fastp.html

    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}${suffix}_fastp.csv
    """
}


process fastp_pre_dehosting {
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_pre-dehosting_fastp.*", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_pre-dehosting_trimmed_R*.fastq.gz", mode: 'copy', enabled: params.publish_trimmed_reads

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_pre-dehosting_fastp.csv")           , emit: metrics
    tuple val(sample_id), path("${sample_id}_pre-dehosting_fastp.json")          , emit: report_json
    tuple val(sample_id), path("${sample_id}_pre-dehosting_fastp.html")          , emit: report_html
    tuple val(sample_id), path("${sample_id}_pre-dehosting_trimmed_R*.fastq.gz") , emit: trimmed_reads

    script:
    worker_threads = task.cpus - 1
    """
    fastp \
      --thread ${worker_threads} \
      -i ${reads[0]} \
      -I ${reads[1]} \
      --cut_tail \
      --trim_poly_g \
      --overrepresentation_analysis \
      --detect_adapter_for_pe \
      -o ${sample_id}_pre-dehosting_trimmed_R1.fastq.gz \
      -O ${sample_id}_pre-dehosting_trimmed_R2.fastq.gz \
      --report_title "fastp report: ${sample_id} (pre-dehosting)" \
      --json ${sample_id}_pre-dehosting_fastp.json \
      --html ${sample_id}_pre-dehosting_fastp.html

    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_pre-dehosting_fastp.json > ${sample_id}_pre-dehosting_fastp.csv
    """
}


process fastp_post_dehosting {
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_post-dehosting_fastp.*", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_post-dehosting_fastp.csv")   , emit: metrics
    tuple val(sample_id), path("${sample_id}_post-dehosting_fastp.json")  , emit: report_json
    tuple val(sample_id), path("${sample_id}_post-dehosting_fastp.html")  , emit: report_html

    script:
    worker_threads = task.cpus - 1
    """
    fastp \
      --thread ${worker_threads} \
      -i ${reads[0]} \
      -I ${reads[1]} \
      --cut_tail \
      --trim_poly_g \
      --overrepresentation_analysis \
      --detect_adapter_for_pe \
      --report_title "fastp report: ${sample_id} (post-dehosting)" \
      --json ${sample_id}_post-dehosting_fastp.json \
      --html ${sample_id}_post-dehosting_fastp.html

    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_post-dehosting_fastp.json > ${sample_id}_post-dehosting_fastp.csv
    """
}
