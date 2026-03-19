#!/usr/bin/env python3
"""
Generate minimal synthetic test data for the PATE pipeline.

Creates:
  - 2 per-locus reference FASTAs (references/locus1.fasta, references/locus2.fasta)
  - Paired-end FASTQ files with embedded heterozygous variants (fastq/)
  - A samplesheet CSV

Usage:
    python generate_test_data.py
"""

import gzip
import os
import random

# ── Configuration ──────────────────────────────────────────────────
SEED = 42
SAMPLE_ID = "sample1"
SPECIES = "TestSpecies"
PLOIDY = 4

NUM_LOCI = 2
LOCUS_LENGTH = 500
READ_LENGTH = 150
FRAGMENT_SIZE = 300
NUM_READ_PAIRS_PER_LOCUS = 120  # ~72x coverage per locus

# SNP positions per locus (well inside the reference, away from edges)
SNP_POSITIONS = {
    "locus1": [100, 200, 300, 400],
    "locus2": [120, 220, 320, 420],
}

# Transition table for generating alt alleles
TRANSITIONS = {"A": "G", "G": "A", "C": "T", "T": "C"}

BASES = ["A", "C", "G", "T"]
QUAL_CHAR = "I"  # Phred 40

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp[b] for b in reversed(seq))


def generate_reference(length, rng):
    """Generate a random DNA sequence of given length."""
    return "".join(rng.choice(BASES) for _ in range(length))


def main():
    rng = random.Random(SEED)

    refs_dir = os.path.join(OUTPUT_DIR, "references")
    fastq_dir = os.path.join(OUTPUT_DIR, "fastq")
    os.makedirs(refs_dir, exist_ok=True)
    os.makedirs(fastq_dir, exist_ok=True)

    # Store reference sequences and alt alleles for read generation
    locus_data = {}

    # ── Generate reference FASTAs ──────────────────────────────────
    for i in range(1, NUM_LOCI + 1):
        locus_name = f"locus{i}"
        ref_seq = list(generate_reference(LOCUS_LENGTH, rng))

        # Determine alt alleles at SNP positions
        snp_pos = SNP_POSITIONS[locus_name]
        alt_alleles = {}
        for pos in snp_pos:
            ref_base = ref_seq[pos]
            alt_alleles[pos] = TRANSITIONS[ref_base]

        locus_data[locus_name] = {
            "ref_seq": "".join(ref_seq),
            "snp_positions": snp_pos,
            "alt_alleles": alt_alleles,
        }

        # Write reference FASTA with sample ID as header
        fasta_path = os.path.join(refs_dir, f"{locus_name}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{SAMPLE_ID}\n{''.join(ref_seq)}\n")
        print(f"  Wrote {fasta_path}")

    # ── Generate paired-end FASTQ reads ────────────────────────────
    r1_path = os.path.join(fastq_dir, f"{SAMPLE_ID}_R1.fastq.gz")
    r2_path = os.path.join(fastq_dir, f"{SAMPLE_ID}_R2.fastq.gz")

    r1_fh = gzip.open(r1_path, "wt")
    r2_fh = gzip.open(r2_path, "wt")

    read_num = 0
    for locus_name, data in locus_data.items():
        ref_seq = data["ref_seq"]
        snp_positions = data["snp_positions"]
        alt_alleles = data["alt_alleles"]

        # Maximum start position so fragment fits within reference
        max_start = LOCUS_LENGTH - FRAGMENT_SIZE
        if max_start < 0:
            max_start = 0

        for pair_idx in range(NUM_READ_PAIRS_PER_LOCUS):
            read_num += 1

            # Random fragment start position
            frag_start = rng.randint(0, max_start)
            frag_end = frag_start + FRAGMENT_SIZE

            # Assign to one of PLOIDY haplotypes (round-robin)
            # Haplotypes 0,1 = reference alleles; 2,3 = alt alleles
            haplotype = pair_idx % PLOIDY

            # Build the fragment sequence with appropriate alleles
            fragment = list(ref_seq[frag_start:frag_end])
            for snp_pos in snp_positions:
                local_pos = snp_pos - frag_start
                if 0 <= local_pos < len(fragment):
                    if haplotype >= PLOIDY // 2:
                        # Alt allele for haplotypes 2,3
                        fragment[local_pos] = alt_alleles[snp_pos]
            fragment = "".join(fragment)

            # R1: forward read from fragment start
            r1_seq = fragment[:READ_LENGTH]
            # R2: reverse complement from fragment end
            r2_seq = reverse_complement(fragment[-READ_LENGTH:])

            qual = QUAL_CHAR * READ_LENGTH
            read_name = f"@{SAMPLE_ID}:{locus_name}:{read_num}"

            r1_fh.write(f"{read_name}/1\n{r1_seq}\n+\n{qual}\n")
            r2_fh.write(f"{read_name}/2\n{r2_seq}\n+\n{qual}\n")

    r1_fh.close()
    r2_fh.close()
    print(f"  Wrote {r1_path} ({read_num} read pairs)")
    print(f"  Wrote {r2_path} ({read_num} read pairs)")

    # ── Generate samplesheet CSV ───────────────────────────────────
    samplesheet_path = os.path.join(OUTPUT_DIR, "samplesheet.csv")
    with open(samplesheet_path, "w") as f:
        f.write("sample,fastq_1,fastq_2,species,ploidy\n")
        f.write(
            f"{SAMPLE_ID},"
            f"fastq/{SAMPLE_ID}_R1.fastq.gz,"
            f"fastq/{SAMPLE_ID}_R2.fastq.gz,"
            f"{SPECIES},"
            f"{PLOIDY}\n"
        )
    print(f"  Wrote {samplesheet_path}")

    # ── Summary ────────────────────────────────────────────────────
    print("\nTest data summary:")
    print(f"  Sample: {SAMPLE_ID}, Species: {SPECIES}, Ploidy: {PLOIDY}")
    print(f"  Loci: {NUM_LOCI} x {LOCUS_LENGTH}bp")
    print(f"  Reads: {read_num} paired-end ({READ_LENGTH}bp), ~{read_num * READ_LENGTH * 2 / (LOCUS_LENGTH * NUM_LOCI):.0f}x coverage")
    for locus_name, data in locus_data.items():
        snps = ", ".join(
            f"{pos}({data['ref_seq'][pos]}->{data['alt_alleles'][pos]})"
            for pos in data["snp_positions"]
        )
        print(f"  {locus_name} SNPs: {snps}")


if __name__ == "__main__":
    main()
