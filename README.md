# Basic Sequence QC

A generic pipeline that can be run on an arbitrary set of Illumina sequence files, regardless of the project or organism of interest.

* Sequence quality information
* Possible contamination

## Analyses

* [`mash`](https://github.com/marbl/Mash): Estimate genome size and depth of coverage
* [`seqtk fqchk`](https://github.com/lh3/seqtk): Measure average sequence quality score, percent of bases above Q30 and GC content.
* [`kraken2`](https://github.com/DerrickWood/kraken2) + [`bracken`](https://github.com/jenniferlu717/Bracken): Taxonomic classification
of reads. Estimation of relative abundances of the most abundant species in each sample.

## Usage

```
nextflow run BCCDC-PHL/basic-sequence-qc \
  [--kraken2_db /path/to/kraken2_db] \
  [--bracken_db /path/to/bracken_db] \
  [--read_length <read_length>] \
  [--seqtk_fqchk_threshold <threshold>] \
  [--mash_sketch_kmer_size <kmer_size>] \
  [--mash_sketch_minimum_copies <copies>] \
  [--prefix 'prefix'] \
  --fastq_input <your illumina run directory> \
  --outdir <output directory>
```

## Output

A single output file in .csv format will be created in the directory specified by `--outdir`. The filename will be `basic_qc_stats.csv`.
If a prefix is provided using the `--prefix` flag, it will be prepended to the output filename, for example: `prefix_basic_qc_stats.csv`.

The output file includes the following headers:

- `sample_id`
- `most_abundant_species_name`
- `most_abundant_species_fraction_total_reads`
- `estimated_genome_size_bp`
- `estimated_depth_coverage`
- `total_bases`
- `average_base_quality`
- `percent_bases_above_q30`
- `percent_gc`
