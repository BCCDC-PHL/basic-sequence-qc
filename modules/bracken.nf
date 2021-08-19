process bracken {

    tag { sample_id }

    errorStrategy 'ignore'

    cpus 1

    input:
      tuple val(sample_id), path(kraken2_report), path(bracken_db)

    output:
      tuple val(sample_id), path("${sample_id}_species_bracken.txt"), path("${sample_id}_species_multiqc_bracken.txt"), path("${sample_id}_species_bracken_abundances.tsv")

    script:
    // MultiQC uses the following regex on the first two lines of a file to identify it as a kraken output:
    // '^\s{1,2}(\d{1,2}\.\d{1,2})\t(\d+)\t(\d+)\t([\dUDKPCOFGS-]{1,3})\t(\d+)\s+(.+)'
    // The output is modified slightly to mimic kraken2 output so that it can be parsed by MultiQC.
    // The original outputs are stored to the output dir, and the modified ones are sent to MultiQC.
    """
    bracken -d ${bracken_db} \
      -i ${kraken2_report} \
      -w ${sample_id}_species_bracken.txt \
      -o ${sample_id}_species_bracken_abundances_unsorted.tsv \
      -r ${params.read_length} \
      -l 'S'
    head -n 1 ${sample_id}_species_bracken_abundances_unsorted.tsv > bracken_abundances_header.tsv
    tail -n+2 ${sample_id}_species_bracken_abundances_unsorted.tsv | sort -t \$'\\t' -nrk 7,7 > ${sample_id}_species_bracken_abundances_data.tsv
    cat bracken_abundances_header.tsv ${sample_id}_species_bracken_abundances_data.tsv > ${sample_id}_species_bracken_abundances.tsv
    sed 's/100\\.00/99\\.99/' ${sample_id}_species_bracken.txt | awk 'NR != 2' | awk '{print " ", \$0}' > ${sample_id}_species_bracken_tmp.txt
    echo -e "  0.01\\t1\\t1\\tU\\t1\\tunclassified" > unclassified_placeholder.tsv
    cat unclassified_placeholder.tsv ${sample_id}_species_bracken_tmp.txt > ${sample_id}_species_multiqc_bracken.txt
    """
}

process abundance_top_n {

    tag { sample_id }

    errorStrategy 'ignore'

    executor 'local'

    cpus 1

    input:
      tuple val(sample_id), path(_), path(_2), path(bracken_abundances)

    output:
      tuple val(sample_id), path("${sample_id}_species_top_*.tsv")

    script:
    def top_n = 5
    """
    bracken_top_n_linelist.py ${bracken_abundances} -n ${top_n} -s ${sample_id} > ${sample_id}_species_top_${top_n}.tsv
    """
}
