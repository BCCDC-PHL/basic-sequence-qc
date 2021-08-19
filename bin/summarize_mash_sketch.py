#!/usr/bin/env python

import argparse
import json

def parse_mash_sketch_output(mash_sketch_output_path):
    output = {}
    with open(mash_sketch_output_path, 'r') as f:
        for line in f:
            if not line.startswith('Writing'):
                fields = [x.strip() for x in line.strip().split(':')]
                key = fields[0].replace(' ', '_').lower()
                value = float(fields[1])
                output[key] = value
                
    return output      


def main(args):
            
    mash_sketch_output = parse_mash_sketch_output(args.mash_sketch_output)

    print(','.join([
        'sample_id',
        'estimated_genome_size_bp',
        'estimated_depth_coverage',
    ]))
    
    print(','.join([
        args.sample_id,
        str(int(mash_sketch_output['estimated_genome_size'])),
        str(round(mash_sketch_output['estimated_coverage'], 3) * 2),
    ]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mash_sketch_output')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
