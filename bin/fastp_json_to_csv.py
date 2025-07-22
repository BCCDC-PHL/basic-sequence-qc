#!/usr/bin/env python

import argparse
import csv
import json
import sys


def seq_is_poly_g(seq):
    """
    Determine if a sequence should be counted as poly-G, based on the following criteria:

    - it includes a contiguous stretch of 10 'G' bases
    - it consists of at least 95% 'G' bases
    """
    if len(seq) < 10:
        return False
    num_g = seq.count('G')
    percent_g = num_g / len(seq)
    if 'GGGGGGGGGG' in seq and percent_g > 0.95:
        return True
    else:
        return False
    

def count_poly_g(fastp_report):
    """
    Count the number of poly-G sequences included in the 'overrepresented_sequences' section of the fastp report.
    """
    poly_g_counts = {}
    r1_overrepresented_seqs_before_filtering = fastp_report['read1_before_filtering']['overrepresented_sequences']

    r1_overrepresented_seqs_after_filtering = fastp_report['read1_after_filtering']['overrepresented_sequences']

    r2_overrepresented_seqs_before_filtering = fastp_report['read2_before_filtering']['overrepresented_sequences']

    r2_overrepresented_seqs_after_filtering = fastp_report['read2_after_filtering']['overrepresented_sequences']

    poly_g_counts['read1_num_poly_g_before_filtering'] = 0
    for seq, count in r1_overrepresented_seqs_before_filtering.items():
        if seq_is_poly_g(seq):
            poly_g_counts['read1_num_poly_g_before_filtering'] += count

    poly_g_counts['read1_num_poly_g_after_filtering'] = 0
    for seq, count in r1_overrepresented_seqs_after_filtering.items():
        if seq_is_poly_g(seq):
            poly_g_counts['read1_num_poly_g_after_filtering'] += count

    poly_g_counts['read2_num_poly_g_before_filtering'] = 0
    for seq, count in r2_overrepresented_seqs_before_filtering.items():
        if seq_is_poly_g(seq):
            poly_g_counts['read2_num_poly_g_before_filtering'] += count

    poly_g_counts['read2_num_poly_g_after_filtering'] = 0
    for seq, count in r2_overrepresented_seqs_after_filtering.items():
        if seq_is_poly_g(seq):
            poly_g_counts['read2_num_poly_g_after_filtering'] += count
    
    return poly_g_counts


def main(args):
    with open(args.fastp_json, 'r') as f:
        fastp_report = json.load(f)

    report_summary = fastp_report['summary']
    before_filtering = report_summary['before_filtering']
    after_filtering = report_summary['after_filtering']
    
    total_reads_before_filtering = before_filtering['total_reads']
    total_reads_after_filtering  = after_filtering['total_reads']
    total_bases_before_filtering = before_filtering['total_bases']
    total_bases_after_filtering  = after_filtering['total_bases']

    read1_mean_length_before_filtering = before_filtering['read1_mean_length']
    read2_mean_length_before_filtering = before_filtering['read2_mean_length']
    read1_mean_length_after_filtering  = after_filtering['read1_mean_length']
    read2_mean_length_after_filtering  = after_filtering['read2_mean_length']

    q20_bases_before_filtering = before_filtering['q20_bases']
    q20_bases_after_filtering  = after_filtering['q20_bases']
    q20_rate_before_filtering  = before_filtering['q20_rate']
    q20_rate_after_filtering   = after_filtering['q20_rate']

    q30_bases_before_filtering = before_filtering['q30_bases']
    q30_bases_after_filtering  = after_filtering['q30_bases']
    q30_rate_before_filtering  = before_filtering['q30_rate']
    q30_rate_after_filtering   = after_filtering['q30_rate']

    gc_content_before_filtering = before_filtering['gc_content']
    gc_content_after_filtering  = after_filtering['gc_content']

    adapter_trimmed_reads = fastp_report['adapter_cutting']['adapter_trimmed_reads']
    adapter_trimmed_bases = fastp_report['adapter_cutting']['adapter_trimmed_bases']

    poly_g_counts = count_poly_g(fastp_report)

    output_fields = [
        'total_reads_before_filtering',
        'total_reads_after_filtering',
        'total_bases_before_filtering',
        'total_bases_after_filtering',
        'read1_mean_length_before_filtering',
        'read1_mean_length_after_filtering',
        'read2_mean_length_before_filtering',
        'read2_mean_length_after_filtering',
        'read1_num_poly_g_before_filtering',
        'read1_num_poly_g_after_filtering',
        'read2_num_poly_g_before_filtering',
        'read2_num_poly_g_after_filtering',
        'q20_bases_before_filtering',
        'q20_bases_after_filtering',
        'q20_rate_before_filtering',
        'q20_rate_after_filtering',
        'q30_bases_before_filtering',
        'q30_bases_after_filtering',
        'q30_rate_before_filtering',
        'q30_rate_after_filtering',
        'gc_content_before_filtering',
        'gc_content_after_filtering',
        'adapter_trimmed_reads',
        'adapter_trimmed_bases',
        'read1_num_poly_g_before_filtering',
        'read1_num_poly_g_after_filtering',
        'read2_num_poly_g_before_filtering',
        'read2_num_poly_g_after_filtering',
    ]

    output_data = {
        'total_reads_before_filtering': total_reads_before_filtering,
        'total_reads_after_filtering': total_reads_after_filtering,
        'total_bases_before_filtering': total_bases_before_filtering,
        'total_bases_after_filtering': total_bases_after_filtering,
        'read1_mean_length_before_filtering': read1_mean_length_before_filtering,
        'read1_mean_length_after_filtering': read1_mean_length_after_filtering,
        'read2_mean_length_before_filtering': read2_mean_length_before_filtering,
        'read2_mean_length_after_filtering': read2_mean_length_after_filtering,
        'q20_bases_before_filtering': q20_bases_before_filtering,
        'q20_bases_after_filtering': q20_bases_after_filtering,
        'q20_rate_before_filtering': q20_rate_before_filtering,
        'q20_rate_after_filtering': q20_rate_after_filtering,
        'q30_bases_before_filtering': q30_bases_before_filtering,
        'q30_bases_after_filtering': q30_bases_after_filtering,
        'q30_rate_before_filtering': q30_rate_before_filtering,
        'q30_rate_after_filtering': q30_rate_after_filtering,
        'gc_content_before_filtering': gc_content_before_filtering,
        'gc_content_after_filtering': gc_content_after_filtering,
        'adapter_trimmed_reads': adapter_trimmed_reads,
        'adapter_trimmed_bases': adapter_trimmed_bases,
    }

    for k, v in poly_g_counts.items():
        output_data[k] = v

    if args.sample_id:
        output_fields = ['sample_id'] + output_fields
        output_data['sample_id'] = args.sample_id

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore', lineterminator='\n')
    writer.writeheader()
    writer.writerow(output_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastp_json')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
