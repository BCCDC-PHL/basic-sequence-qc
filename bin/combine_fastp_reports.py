#!/usr/bin/env python3

import argparse
import csv
import json
import sys

from pathlib import Path

def parse_fastp_csv(fastp_csv_path: Path):
    """
    """
    fastp_report = {}
    int_fields = [
        'total_reads_before_filtering',
        'total_reads_after_filtering',
        'total_bases_before_filtering',
        'total_bases_after_filtering',
        'read1_mean_length_before_filtering',
        'read1_mean_length_after_filtering',
        'read2_mean_length_before_filtering',
        'read2_mean_length_after_filtering',
        'q20_bases_before_filtering',
        'q30_bases_before_filtering',
        'q20_bases_after_filtering',
        'q30_bases_after_filtering',
        'adapter_trimmed_reads',
        'adapter_trimmed_bases',
        'read1_num_poly_g_before_filtering',
        'read1_num_poly_g_after_filtering',
        'read2_num_poly_g_before_filtering',
        'read2_num_poly_g_after_filtering',        
    ]
    float_fields = [
        'q20_rate_before_filtering',
        'q30_rate_before_filtering',
        'q20_rate_after_filtering',
        'q30_rate_after_filtering',
        'gc_content_before_filtering',
        'gc_content_after_filtering',
    ]
    with open(fastp_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for field, value in row.items():
                if field in int_fields:
                    try:
                        fastp_report[field] = int(value)
                    except ValueError as e:
                        fastp_report[field] = None
                elif field in float_fields:
                    try:
                        fastp_report[field] = float(value)
                    except ValueError as e:
                        fastp_report[field] = None
                else:
                    fastp_report[field] = value

    return fastp_report


def combine_fastp_reports(pre_dehosting_report, post_dehosting_report):
    """
    """
    combined_report = {}

    combined_report['sample_id'] = pre_dehosting_report['sample_id']

    for key, value in pre_dehosting_report.items():
        if key.endswith('before_filtering'):
            combined_report[key] = value
    for key, value in post_dehosting_report.items():
        if key.endswith('after_filtering'):
            combined_report[key] = value

    pre_dehosting_adapter_trimmed_reads = pre_dehosting_report.get('adapter_trimmed_reads', 0)
    post_dehosting_adapter_trimmed_reads = post_dehosting_report.get('adapter_trimmed_reads', 0)
    combined_report['adapter_trimmed_reads'] = pre_dehosting_adapter_trimmed_reads + post_dehosting_adapter_trimmed_reads

    pre_dehosting_adapter_trimmed_bases = pre_dehosting_report.get('adapter_trimmed_bases', 0)
    post_dehosting_adapter_trimmed_bases = post_dehosting_report.get('adapter_trimmed_bases', 0)
    combined_report['adapter_trimmed_bases'] = pre_dehosting_adapter_trimmed_bases + post_dehosting_adapter_trimmed_bases
        
    return combined_report


def main(args):
    pre_dehosting_report = parse_fastp_csv(Path(args.pre_dehosting))
    post_dehosting_report = parse_fastp_csv(Path(args.post_dehosting))

    combined_report = combine_fastp_reports(pre_dehosting_report, post_dehosting_report)

    output_fieldnames = [
        'sample_id',
        'total_reads_before_filtering',
        'total_reads_after_filtering',
        'total_bases_before_filtering',
        'total_bases_after_filtering',
        'read1_mean_length_before_filtering',
        'read1_mean_length_after_filtering',
        'read2_mean_length_before_filtering',
        'read2_mean_length_after_filtering',
        'q20_bases_before_filtering',
        'q30_bases_after_filtering',
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
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()
    writer.writerow(combined_report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pre-dehosting', required=True, help='Path to pre-dehosting fastp report')
    parser.add_argument('--post-dehosting', required=True, help='Path to post-dehosting fastp report')
    args = parser.parse_args()
    main(args)
