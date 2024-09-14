# Basic Sequence QC

A generic pipeline that can be run on an arbitrary set of Illumina sequence files, regardless of the project or organism of interest.

* Sequence quality information

## Analyses

* [`fastp`](https://github.com/OpenGene/fastp): Collect sequence QC stats

## Usage

```
nextflow run BCCDC-PHL/basic-sequence-qc \
  [--prefix 'prefix'] \
  --fastq_input <your fastq input directory> \
  --outdir <output directory>
```

Alternatively, a `sample_sheet.csv` file can be provided, with fields:

```
ID
R1
R2
```

For example:
```csv
ID,R1,R2
SAMPLE-01,/path/to/SAMPLE-01_R1.fastq.gz,/path/to/SAMPLE-01_R2.fastq.gz
SAMPLE-02,/path/to/SAMPLE-02_R1.fastq.gz,/path/to/SAMPLE-02_R2.fastq.gz
SAMPLE-03,/path/to/SAMPLE-03_R1.fastq.gz,/path/to/SAMPLE-03_R2.fastq.gz
```

The `sample_sheet.csv` file can be provided using the `--sample_sheet_input` flag as follows:

```
nextflow run BCCDC-PHL/basic-sequence-qc \
  [--prefix 'prefix'] \
  --sample_sheet_input <your sample_sheet.csv file> \
  --outdir <output directory>
```

## Dehosting

This pipeline supports an optional dehosting...


## Output

A single output file in .csv format will be created in the directory specified by `--outdir`. The filename will be `basic_qc_stats.csv`.
If a prefix is provided using the `--prefix` flag, it will be prepended to the output filename, for example: `prefix_basic_qc_stats.csv`.

The output file includes the following headers:

```
sample_id
total_reads_before_filtering
total_reads_after_filtering
total_bases_before_filtering
total_bases_after_filtering
read1_mean_length_before_filtering
read1_mean_length_after_filtering
read2_mean_length_before_filtering
read2_mean_length_after_filtering
q20_bases_before_filtering
q20_bases_after_filtering
q20_rate_before_filtering
q20_rate_after_filtering
q30_bases_before_filtering
q30_bases_after_filtering
q30_rate_before_filtering
q30_rate_after_filtering
gc_content_before_filtering
gc_content_after_filtering
adapter_trimmed_reads
adapter_trimmed_bases
```
