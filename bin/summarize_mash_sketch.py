#!/usr/bin/env python

import argparse
import csv
import json
import sys

def parse_mash_sketch_output(mash_sketch_output_path):
    output = {}
    with open(mash_sketch_output_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('Writing') and not line.startswith('ERROR') and line != "":
                fields = [x.strip() for x in line.split(':')]
                key = fields[0].replace(' ', '_').lower()
                value = float(fields[1])
                output[key] = value
                
    return output      


def main(args):
            
    mash_sketch_output = parse_mash_sketch_output(args.mash_sketch_output)

    mash_sketch_output['sample_id'] = args.sample_id

    if 'estimated_coverage' in mash_sketch_output:
        mash_sketch_output['estimated_coverage'] = round(mash_sketch_output['estimated_coverage'], 3) * 2

    fieldname_translation = {
        'sample_id': 'sample_id',
        'estimated_coverage': 'estimated_depth_coverage',
        'estimated_genome_size': 'estimated_genome_size_bp',
    }

    translated_mash_sketch_output = {}    
    for k, v in fieldname_translation.items():
        if k in mash_sketch_output:
            translated_mash_sketch_output[v] = mash_sketch_output[k]
        
    output_fieldnames = [
        'sample_id',
        'estimated_genome_size_bp',
        'estimated_depth_coverage',
    ]

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    writer.writerow(translated_mash_sketch_output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mash_sketch_output')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
