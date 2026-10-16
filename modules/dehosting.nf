process dehost {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_dehosted_R{1,2}.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_hostile.log.json", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(hostile_cache_dir)

    output:
    tuple val(sample_id), path("${sample_id}_hostile.log.json"), emit: hostile_log
    tuple val(sample_id), path("${sample_id}_dehosted*.fastq.gz"), emit: dehosted_reads
    tuple val(sample_id), path("${sample_id}_dehost_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: dehost\\n"                                            >> ${sample_id}_dehost_provenance.yml
    printf -- "  tools:\\n"                                                          >> ${sample_id}_dehost_provenance.yml
    printf -- "    - tool_name: hostile\\n"                                          >> ${sample_id}_dehost_provenance.yml
    printf -- "      tool_version: \$(hostile --version 2>&1 | tail -n 1)\\n"         >> ${sample_id}_dehost_provenance.yml
    printf -- "      parameters:\\n"                                                 >> ${sample_id}_dehost_provenance.yml
    printf -- "        - parameter: --index\\n"                                      >> ${sample_id}_dehost_provenance.yml
    printf -- "          value: ${params.dehosting_index}\\n"                        >> ${sample_id}_dehost_provenance.yml
    printf -- "        - parameter: --threads\\n"                                    >> ${sample_id}_dehost_provenance.yml
    printf -- "          value: ${task.cpus}\\n"                                     >> ${sample_id}_dehost_provenance.yml

    export HOSTILE_CACHE_DIR=${hostile_cache_dir}

    hostile clean \
	--threads ${task.cpus} \
	--fastq1 ${reads[0]} \
	--fastq2 ${reads[1]} \
	--index ${params.dehosting_index} \
	--output . \
	> ${sample_id}_hostile.log.json

    mv ${sample_id}*.clean_1.fastq.gz ${sample_id}_dehosted_R1.fastq.gz
    mv ${sample_id}*.clean_2.fastq.gz ${sample_id}_dehosted_R2.fastq.gz
    """
}


process combine_fastp_reports {

    tag { sample_id }

    input:
    tuple val(sample_id), path(fastp_pre_dehosting), path(fastp_post_dehosting)

    output:
    tuple val(sample_id), path("${sample_id}_fastp_combined.csv"), emit: metrics

    script:
    """
    combine_fastp_reports.py \
	--pre-dehosting ${fastp_pre_dehosting} \
	--post-dehosting ${fastp_post_dehosting} \
	> ${sample_id}_fastp_combined.csv
    """
}
