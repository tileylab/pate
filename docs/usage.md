# tileylab/pate: Usage

## Introduction

tileylab/pate is a Nextflow pipeline for recovering phased haplotype sequences from
short-read target enrichment data, particularly for polyploid organisms. It reimplements
the [PATE (Phased Alleles from Target Enrichment)](https://github.com/gtiley/Phasing)
Perl pipeline as a portable, containerized Nextflow workflow.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like
to analyse before running the pipeline. Use the `--input` parameter to specify its location.

```csv
sample,fastq_1,fastq_2,ploidy,species
KJM225,/path/to/KJM225.R1.fq.gz,/path/to/KJM225.R2.fq.gz,2,Cyperus
KJM226,/path/to/KJM226.R1.fq.gz,/path/to/KJM226.R2.fq.gz,4,Cyperus
KJM300,/path/to/KJM300.R1.fq.gz,/path/to/KJM300.R2.fq.gz,6,Carex
```

| Column    | Description                                                      |
| --------- | ---------------------------------------------------------------- |
| `sample`  | Unique sample identifier. Must match headers in reference FASTAs |
| `fastq_1` | Path to read 1 FASTQ file (gzipped)                             |
| `fastq_2` | Path to read 2 FASTQ file (gzipped)                             |
| `ploidy`  | Integer ploidy level for this sample (e.g., 2, 4, 6)            |
| `species` | Species or group identifier for matching reference sequences     |

## Reference files

The `--references` parameter should point to a directory containing per-locus FASTA files.
Each file contains one sequence per sample/species:

```
references/
  locus001.fasta
  locus002.fasta
  ...
```

Each FASTA file should contain sequences with headers matching either the `sample` or
`species` column from your samplesheet:

```
>KJM225
ACGTACGTACGT...
>KJM226
ACGTACGTACGT...
```

These are typically supercontig assemblies from tools like HybPiper.

## Running the pipeline

The typical command for running the pipeline is:

```bash
nextflow run tileylab/pate \
    -profile docker \
    --input samplesheet.csv \
    --references /path/to/references \
    --outdir results
```

### Key parameters

| Parameter                  | Default     | Description                                      |
| -------------------------- | ----------- | ------------------------------------------------ |
| `--input`                  | (required)  | Path to samplesheet CSV                          |
| `--references`             | (required)  | Path to directory of per-locus reference FASTAs   |
| `--outdir`                 | (required)  | Output directory                                 |
| `--remove_duplicates`      | `false`     | Mark PCR duplicates with Picard MarkDuplicates   |
| `--unique_only`            | `true`      | Output only unique haplotype sequences           |
| `--output_expected_dosage` | `false`     | Output ploidy-count copies for invariant loci    |
| `--genotype_mode`          | `consensus` | `consensus` or `iupac` for unphased variants     |

### GATK variant filter parameters

| Parameter               | Default | Description                 |
| ----------------------- | ------- | --------------------------- |
| `--variant_filter_qd`   | `2.0`   | Quality by depth threshold  |
| `--variant_filter_fs`   | `60.0`  | Fisher strand bias          |
| `--variant_filter_mq`   | `40.0`  | Mapping quality             |
| `--variant_filter_rprs` | `-8.0`  | Read position rank sum      |
| `--variant_filter_af_lo`| `0.025` | Allele frequency lower      |
| `--variant_filter_af_hi`| `0.975` | Allele frequency upper      |
| `--variant_filter_dp`   | `10`    | Minimum depth               |

## Profiles

Use `-profile` to select a container engine and optional configurations:

- `docker` - Docker containers
- `singularity` - Singularity containers
- `podman` - Podman containers
- `conda` - Conda environments
- `test` - Minimal test dataset

Multiple profiles can be combined: `-profile test,docker`
