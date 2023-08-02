#!/usr/bin/env python

import argparse
import csv
import json
import sys

def parse_seqtk_fqchk_output(seqtk_fqchk_output_path, quality_threshold):
    output = []
    with open(seqtk_fqchk_output_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            parsed_row = {}
            if row['position'] == 'ALL':
                percent_above_header = 'percent_bases_above_q' + str(quality_threshold)
                parsed_row['num_bases'] = int(row['num_bases'])
                parsed_row['average_q'] = float(row['average_q'])
                parsed_row[percent_above_header] = float(row[percent_above_header])
                percent_g = float(row['percent_g'])
                percent_c = float(row['percent_c'])
                parsed_row['percent_gc'] = percent_g + percent_c
                output.append(parsed_row)

    return output      


def main(args):
    header = []
    with open(args.seqtk_fqchk_output,'r') as f:
        header = f.readline().strip().split(',')

    quality_threshold = 0
    for h in header:
        if h.startswith('percent_bases_above_q'):
            quality_threshold = int(h.split('percent_bases_above_q')[1])
            
    seqtk_fqchk_output = parse_seqtk_fqchk_output(args.seqtk_fqchk_output, quality_threshold)
    
    total_bases = sum([x['num_bases'] for x in seqtk_fqchk_output])

    percent_above_header = 'percent_bases_above_q' + str(quality_threshold)

    if total_bases == 0:
        overall_percent_gc = 0.0
        overall_average_q = 0.0
        overall_percent_above_threshold = 0.0
    else:
        overall_percent_gc = sum([x['percent_gc'] * x['num_bases'] for x in seqtk_fqchk_output]) / total_bases
        overall_average_q = sum([x['average_q'] * x['num_bases'] for x in seqtk_fqchk_output]) / total_bases
        overall_percent_above_threshold = sum([x[percent_above_header] * x['num_bases'] for x in seqtk_fqchk_output]) / total_bases    

    print(','.join([
        'sample_id',
        'percent_gc',
        'total_bases',
        'average_base_quality',
        percent_above_header,
    ]))
    
    print(','.join([
        args.sample_id,
        str(round(overall_percent_gc, 3)),
        str(total_bases),
        str(round(overall_average_q, 3)),
        str(round(overall_percent_above_threshold, 3)),
    ]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('seqtk_fqchk_output')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
