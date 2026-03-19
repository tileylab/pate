# tileylab/pate

## Phased Alleles from Target Enrichment (PATE)

A Nextflow DSL2 pipeline for recovering phased haplotype sequences from short-read target enrichment data of polyploid organisms. PATE maps reads to per-sample reference sequences, calls variants with GATK4 HaplotypeCaller at user-specified ploidy, phases haplotypes with H-PoPG, and produces per-locus phased allele sequences.

This pipeline is a Nextflow reimplementation of the [original Perl-based PATE pipeline](https://github.com/gtiley/Phasing), following [nf-core](https://nf-co.re) conventions for reproducibility, containerization, and portability.

## Pipeline summary

1. Read quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Per-sample reference preparation from per-locus FASTA files
3. Read mapping ([BWA-MEM](https://github.com/lh3/bwa)) with BAM merging via GATK MergeBamAlignment
4. Optional PCR duplicate removal ([Picard MarkDuplicates](https://broadinstitute.github.io/picard/))
5. BAM statistics ([samtools](https://www.htslib.org/))
6. Variant calling at user-specified ploidy ([GATK4 HaplotypeCaller](https://gatk.broadinstitute.org/))
7. Variant filtering: SNP selection, allele fraction annotation, hard filtering, biallelic filtering
8. IUPAC consensus generation ([GATK4 FastaAlternateReferenceMaker](https://gatk.broadinstitute.org/))
9. Per-locus BAM splitting and haplotype phasing ([H-PoPG](https://github.com/MinzhuXie/H-PoPG))
10. Phased consensus sequence construction and summary statistics
11. Aggregate report ([MultiQC](http://multiqc.info/))

## Quick start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) (`>=25.04.0`)

2. Install [Docker](https://docs.docker.com/engine/installation/), [Singularity](https://www.sylabs.io/guides/3.0/user-guide/), or [Conda](https://conda.io/miniconda.html)

3. Run the pipeline:

   ```bash
   nextflow run tileylab/pate \
       -profile docker \
       --input samplesheet.csv \
       --references /path/to/references \
       --outdir results
   ```

   On Apple Silicon (M1/M2/M3/M4), add the `emulate_amd64` profile:

   ```bash
   nextflow run tileylab/pate \
       -profile docker,emulate_amd64 \
       --input samplesheet.csv \
       --references /path/to/references \
       --outdir results
   ```

## Samplesheet

The samplesheet is a CSV file with information about the samples to process. It must have the following columns:

| Column    | Description                                                                 |
|-----------|-----------------------------------------------------------------------------|
| `sample`  | Sample identifier (no spaces). Used as FASTA header match and file prefix.  |
| `fastq_1` | Path to read 1 FASTQ file (gzipped).                                       |
| `fastq_2` | Path to read 2 FASTQ file (gzipped).                                       |
| `species` | Species or population identifier. Used as a fallback match in per-locus reference FASTAs. |
| `ploidy`  | Integer ploidy level (e.g. `2` for diploid, `4` for tetraploid).            |

Example:

```csv
sample,fastq_1,fastq_2,species,ploidy
IndividualA,/data/IndividualA_R1.fastq.gz,/data/IndividualA_R2.fastq.gz,SpeciesX,4
IndividualB,/data/IndividualB_R1.fastq.gz,/data/IndividualB_R2.fastq.gz,SpeciesY,2
```

## Reference FASTAs

The `--references` directory should contain one FASTA file per locus. Each file contains one sequence per sample (or species), with the FASTA header matching the `sample` or `species` column from the samplesheet.

```
references/
  locus1.fasta    # >IndividualA\nACGT...\n>IndividualB\nACGT...
  locus2.fasta
  locus3.fasta
```

Rules for reference files:
- One FASTA file per locus, named `<locus_name>.fasta` (or `.fa`, `.fna`)
- The locus filename (without extension) becomes the contig name in the combined reference
- FASTA headers should match sample IDs; species names are used as a fallback
- No spaces in locus names or FASTA headers
- Sequences should not contain line breaks within the sequence data

## Pipeline parameters

### Input/output options

| Parameter      | Default  | Description                                                    |
|----------------|----------|----------------------------------------------------------------|
| `--input`      | *required* | Path to samplesheet CSV.                                     |
| `--references` | *required* | Path to directory containing per-locus reference FASTA files. |
| `--outdir`     | *required* | Output directory for results.                                 |
| `--email`      | `null`   | Email address for completion summary.                          |

### Processing options

| Parameter                  | Default     | Description |
|----------------------------|-------------|-------------|
| `--remove_duplicates`      | `false`     | Mark and remove PCR duplicates with Picard MarkDuplicates. Leave `false` for target enrichment data; set `true` for non-enriched (e.g. whole-genome) libraries. |
| `--unique_only`            | `true`      | Output only unique phased haplotype sequences. When `false`, the expected number of alleles (matching ploidy) are output even if some are identical. |
| `--output_expected_dosage` | `false`     | Output ploidy-count copies of the reference sequence for invariant loci (loci with no variants). May be useful when calculating allele frequencies. |
| `--genotype_mode`          | `consensus` | How to handle unphased variants. `consensus`: replace unphased positions with `N`. `iupac`: replace unphased positions with IUPAC ambiguity codes. |

When multiple haplotype blocks are recovered for a locus, the pipeline retains phasing from the longest block only. Variants in shorter blocks are treated according to `--genotype_mode`.

### Variant filtering options

These thresholds control GATK VariantFiltration hard filters applied after variant calling. Variants failing any filter are excluded from phasing.

| Parameter               | Default | Filter expression           | Description |
|-------------------------|---------|-----------------------------|-------------|
| `--variant_filter_qd`   | `2.0`   | `QD < 2.0`                 | Quality by depth. Low values indicate low-confidence variants. |
| `--variant_filter_fs`   | `60.0`  | `FS > 60.0`                | Fisher strand bias. High values indicate strand bias artifacts. |
| `--variant_filter_mq`   | `40.0`  | `MQ < 40.0`                | Root mean square mapping quality. Low values indicate poor mapping. |
| `--variant_filter_rprs` | `-8.0`  | `ReadPosRankSum < -8.0`    | Read position rank sum. Extreme negative values indicate positional bias. |
| `--variant_filter_af_lo`| `0.025` | `AF < 0.025`               | Allele frequency lower bound. Removes very rare alleles likely from error. |
| `--variant_filter_af_hi`| `0.975` | `AF > 0.975`               | Allele frequency upper bound. Removes near-fixed alleles (likely reference error). |
| `--variant_filter_dp`   | `10`    | `DP < 10`                  | Minimum read depth. Variants with fewer supporting reads are filtered. |

## Output directories

```
results/
  fastqc/           # FastQC reports per sample
  references/       # Per-sample concatenated reference FASTAs
  bwa/              # BWA index files and aligned BAMs
  picard/           # Merged BAMs and unmapped BAMs
  samtools/         # BAM index and alignment statistics
  bamtools/         # Per-locus split BAMs
  gatk4/            # HaplotypeCaller VCFs, annotated and filtered VCFs, sequence dictionaries
  select/           # SNP-selected and final filtered VCFs
  iupac/            # IUPAC consensus FASTA (heterozygous positions as ambiguity codes)
  phasing/          # H-PoPG phase output and log files per locus
  PHASED/           # Phased haplotype FASTA sequences per locus
  GENOTYPE/         # Unphased genotype sequences (consensus or IUPAC-coded)
  PICKONE/          # One randomly chosen haplotype per locus per sample
  summary_stats/    # Per-sample phasing summaries and averagePhasingStats.txt
  multiqc/          # Aggregate QC report
  pipeline_info/    # Execution reports, timelines, traces, and DAGs
```

### Key output files

- **`PHASED/*.phased.fasta`** — Per-locus FASTA files containing all phased haplotype sequences. Each sequence is named `>sampleID__locusName__alleleNumber`.
- **`GENOTYPE/*.genotype.fasta`** — Per-locus FASTA with the unphased genotype sequence, where heterozygous positions use consensus (`N`) or IUPAC codes depending on `--genotype_mode`.
- **`PICKONE/*.pickone.fasta`** — Per-locus FASTA with a single randomly chosen allele per sample. Useful for phylogenetic network analyses where only one sequence per individual is desired.
- **`summary_stats/averagePhasingStats.txt`** — Tab-delimited summary across all samples and loci.

### Summary statistics columns

| Column         | Description |
|----------------|-------------|
| `INDIVIDUAL`   | Sample ID |
| `NLOCI`        | Number of loci assembled |
| `NVARLOCI`     | Number of loci with at least one variant |
| `NINVLOCI`     | Number of invariant loci |
| `NPHASELOCI`   | Number of loci that were phased |
| `AVG_LENGTH`   | Average locus length (bp) |
| `AVG_NVAR`     | Average number of variants per locus (including invariant loci) |
| `AVG_HET`      | Average per-base heterozygosity across all loci |
| `AVG_NBLOCKS`  | Average number of phasing blocks per phased locus |
| `AVG_LONGESTBL`| Average number of variants in the longest phasing block |
| `1, 2, 3, ...` | Number of loci with that many phased alleles (columns extend to the maximum ploidy) |

## Profiles

| Profile          | Description |
|------------------|-------------|
| `docker`         | Run with Docker containers |
| `singularity`    | Run with Singularity containers |
| `conda`          | Run with Conda environments |
| `emulate_amd64`  | Add `--platform=linux/amd64` to Docker. Required for Apple Silicon Macs since bioinformatics containers are typically amd64-only. |
| `test`           | Run with bundled minimal test dataset |

## Software dependencies

All dependencies are provided automatically through containers or Conda. The pipeline uses:

- [BWA](https://github.com/lh3/bwa) (Li & Durbin 2009)
- [samtools/htslib](https://www.htslib.org/) (Li et al. 2009)
- [GATK4](https://gatk.broadinstitute.org/) (McKenna et al. 2010)
- [H-PoPG](https://github.com/MinzhuXie/H-PoPG) (Xie et al. 2016)
- [Picard](https://broadinstitute.github.io/picard/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](http://multiqc.info/)

## Citations

If you use this pipeline, please cite:

- Tiley GP, Crowl AA, Manos PS, Sessa EB, Solis-Lemus C, Yoder AD, Burleigh JG. 2021. Phasing alleles improves network inference with allopolyploids. *bioRxiv*. doi: [10.1101/2021.05.04.442457](https://doi.org/10.1101/2021.05.04.442457)
- Crowl AA, Fritsch PW, Tiley GP, Lynch NP, Ranney TG, Ashrafi H, Manos PS. 2022. A first complete phylogenomic hypothesis for diploid blueberries (*Vaccinium* section *Cyanococcus*). *American Journal of Botany*.
- Tiley GP, Crowl AA, Almary TOM, Luke WRQ, Solofondranohatra CL, Besnard G, Lehmann CER, Yoder AD, Vorontsova MS. 2023. Genetic variation in *Loudetia simplex* supports the presence of ancient grasslands in Madagascar. *bioRxiv*. doi: [10.1101/2023.04.07.536094](https://doi.org/10.1101/2023.04.07.536094)

Please also cite the tools used in the pipeline:

- Li H, Durbin R. 2009. Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics* 25:1754-1760.
- Li H, Handsaker B, Wysoker A, et al. 2009. The sequence alignment/map format and SAMtools. *Bioinformatics* 25:2078-2079.
- McKenna A, Hanna M, Banks E, et al. 2010. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res* 20:1297-1303.
- Xie M, Wu Q, Wang J, Jiang T. 2016. H-PoP and H-PoPG: heuristic partitioning algorithms for single individual haplotyping of polyploids. *Bioinformatics* 32:3735-3744.

## Credits

tileylab/pate was originally written by [George P. Tiley](https://github.com/gtiley).

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community.
