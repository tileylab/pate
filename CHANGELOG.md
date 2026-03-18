# tileylab/pate: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of tileylab/pate - a Nextflow implementation of the PATE
(Phased Alleles from Target Enrichment) pipeline.

### `Added`

- Species mode: per-sample reference preparation, read mapping, variant calling,
  GATK filtering, H-PoPG haplotype phasing, and allele extraction
- Samplesheet-based input with per-sample ploidy and species metadata
- Per-locus parallelization for phasing and allele extraction
- Output: PHASED, GENOTYPE, and PICKONE FASTA files per locus
- Summary phasing statistics per individual and global averages
- MultiQC report integration (FastQC, samtools stats)
- Container support: Docker, Singularity, Podman, Apptainer
- Configurable GATK variant filter parameters
- Optional duplicate marking via Picard MarkDuplicates
