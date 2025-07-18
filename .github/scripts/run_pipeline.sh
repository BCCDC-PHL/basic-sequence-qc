#!/bin/bash

set -eo pipefail

nextflow run main.nf \
	 -profile "${PROFILE}" \
	 --fastq_input .github/data/fastq \
	 --outdir .github/data/test_output \
	 --prefix test \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv
