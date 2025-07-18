#!/usr/bin/env python3

import argparse
import csv
import glob
import json
import os
import urllib.request

from jsonschema import validate
import yaml


def check_provenance_format_valid(provenance_files, schema):
    """
    Check that the provenance files are valid according to the schema.
    """
    for provenance_file in provenance_files:
        with open(provenance_file) as f:
            try:
                provenance = yaml.load(f, Loader=yaml.BaseLoader)
                validate(provenance, schema)
            except Exception as e:
                print(f"Error validating {provenance_file}: {e}")
                exit(1)
                return False

    return True

def check_expected_files_exist(output_dir, prefix="test"):
    """
    Check that the expected files exist in the output directory.

    :param output_dir: Path to the output directory
    :param sample_ids: List of sample IDs
    :return: True if all expected files exist, False otherwise
    :rtype: bool
    """
    expected_files = [
        f"{prefix}_basic_qc_stats.csv",
    ]

    for expected_file in expected_files:
        expected_file_path = os.path.join(output_dir, expected_file)
        if not os.path.exists(expected_file_path):
            print(f"Expected file {expected_file_path} not found")
            return False

    return True


def main(args):

    output_dir = os.path.dirname(args.output)
    os.makedirs(output_dir, exist_ok=True)

    # TODO: Add more tests
    tests = [
        {
            "test_name": "all_expected_files_exist",
            "test_passed": check_expected_files_exist(args.pipeline_outdir, args.prefix),
        },
    ]

    output_fields = [
        "test_name",
        "test_result"
    ]

    output_path = args.output
    with open(output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=output_fields, extrasaction='ignore')
        writer.writeheader()
        for test in tests:
            if test["test_passed"]:
                test["test_result"] = "PASS"
            else:
                test["test_result"] = "FAIL"
            writer.writerow(test)

    for test in tests:
        if not test['test_passed']:
            exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check outputs')
    parser.add_argument('--pipeline-outdir', type=str, help='Path to the pipeline output directory')
    parser.add_argument('--prefix', type=str, default='test', help='Prefix used for test pipeline outputs')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    main(args)
