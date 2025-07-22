# Basic Sequence QC

A generic pipeline that can be run on an arbitrary set of Illumina sequence files, regardless of the project or organism of interest.

* Sequence quality information

## Analyses

* [`fastp`](https://github.com/OpenGene/fastp): Collect sequence QC stats
* [`hostile`](https://github.com/bede/hostile): Removal of host (human) reads (optional)

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

### Publishing Trimmed Reads

By default, the trimmed reads are not 'published' to the output directory. To enable publishing the trimmed reads,
add the `--publish_trimmed_reads` flag:

```
nextflow run BCCDC-PHL/basic-sequence-qc \
  [--prefix 'prefix'] \
  --publish_trimmed_reads \
  --sample_sheet_input <your sample_sheet.csv file> \
  --outdir <output directory>
```

## Dehosting

This pipeline supports an optional dehosting step that can be used to remove human-derived sequence reads.
Removal of host reads from other organisms is not currently supported.

### Setup

Before using this pipeline to perform dehosting, the human reference genome(s) should be downloaded.
Dehosting is performed by attempting to align reads against the human genome and removing reads
that align with sufficient quality and specificity.

1. Activate the conda environment that is used by this pipeline for dehosting.

```
conda activate basic-sequence-qc-dehosting-eeb6f6f753e1d37317d1d06b4e3141d5
```

2. Download (fetch) the default reference genome (`human-t2t-hla`)

```
hostile fetch
```

The reference will be downloaded to:

```
~/.local/share/hostile
```

If alternative versions of the reference genome are needed, download them 

```
hostile fetch --name <REF_NAME>
```

Details on available references can be found on [the bede/hostile README](https://github.com/bede/hostile?tab=readme-ov-file#indexes).
The list of available reference names can be found by running:

```
hostile fetch --list
```

...which should return a list like this:

```
human-t2t-hla
human-t2t-hla-argos985
human-t2t-hla-argos985-mycob140
human-t2t-hla.rs-viral-202401_ml-phage-202401
human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401
```

3. Deactivate the conda env.

```
conda deactivate
```

### Performing Dehosting

To include the dehosting step in the pipeline analysis, add the `--dehost` flag:

```
nextflow run BCCDC-PHL/basic-sequence-qc \
  -profile conda \
  --cache ~/.conda/envs \
  --fastq_input /path/to/fastq_input \
  --dehost \
  --outdir /path/to/output-dir
```

By default, the `human-t2t-hla` reference will be used. To use an alternative reference, include the `--dehosting_index` flag along
with the reference name:

```
nextflow run BCCDC-PHL/basic-sequence-qc \
  -profile conda \
  --cache ~/.conda/envs \
  --fastq_input /path/to/fastq_input \
  --dehost \
  --dehosting_index human-t2t-hla-argos985-mycob140 \
  --outdir /path/to/output-dir
```

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

### Dehosted Reads

If dehosting is performed, dehosted reads will be deposited under the directory supplied for the `--outdir` param, in a sub-directory
for each sample, named using the sample ID. For example:

```
outdir
├── sample-01
│   ├── sample-01_dehosted_R1.fastq.gz
│   ├── sample-01_dehosted_R2.fastq.gz
│   ├── sample-01_dehosted_fastp.csv
│   ├── sample-01_fastp.csv
│   └── sample-01_hostile.log.json
├── sample-02
│   ├── sample-02_dehosted_R1.fastq.gz
│   ├── sample-02_dehosted_R2.fastq.gz
│   ├── sample-02_dehosted_fastp.csv
│   ├── sample-02_fastp.csv
│   └── sample-02_hostile.log.json
├── sample-03
│   ├── ...
...
```
