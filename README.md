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
