#!/usr/bin/env python3
"""
Generate minimal synthetic test data for the PATE pipeline.

Produces a small but *representative* dataset so the end-to-end outputs can be
validated and regressions (e.g. the cross-locus consensus bug) are exercised:

  - Multiple samples with mixed ploidy and alt-allele dosage
  - Several per-locus reference FASTAs, each carrying one record per sample
    (header = sample id), all sharing the locus reference sequence
  - Paired-end FASTQ per sample with embedded heterozygous variants
  - A samplesheet CSV (absolute FASTQ paths, so nf-schema resolves them
    regardless of launch directory)

Deterministic via SEED. Usage:
    python generate_test_data.py
"""

import gzip
import os
import random

# ── Configuration ──────────────────────────────────────────────────
SEED = 42
SPECIES = "TestSpecies"

# sample id -> ploidy + number of haplotypes carrying the ALT allele at each SNP
# (alt_dosage < ploidy keeps every SNP heterozygous)
SAMPLES = [
    {"id": "sample1", "ploidy": 4, "alt_dosage": 2},  # tetraploid, balanced 2:2
    {"id": "sample2", "ploidy": 2, "alt_dosage": 1},  # diploid, 1:1
    {"id": "sample3", "ploidy": 4, "alt_dosage": 3},  # tetraploid, skewed 3:1
]

# locus name -> reference length (varied, to exercise the LENGTH column)
LOCI = {
    "locus1": 500,
    "locus2": 400,
    "locus3": 600,
    "locus4": 450,
    "locus5": 550,
}

READ_LENGTH = 150
FRAGMENT_SIZE = 280               # fits within the shortest locus (400 bp)
NUM_READ_PAIRS_PER_LOCUS = 100    # per sample per locus (~60x)
NUM_SNPS = 6                      # SNPs per locus

# Transition table for generating alt alleles
TRANSITIONS = {"A": "G", "G": "A", "C": "T", "T": "C"}
BASES = ["A", "C", "G", "T"]
QUAL_CHAR = "I"  # Phred 40

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp[b] for b in reversed(seq))


def snp_positions(length, locus_index):
    """Deterministic SNP positions: clustered within a fragment (so they phase)
    but shifted per locus so positions differ across loci."""
    offset = 40 + 15 * locus_index
    positions = [offset + 45 * i for i in range(NUM_SNPS)]
    # keep away from the locus edges
    return [p for p in positions if 20 <= p < length - 20]


def main():
    rng = random.Random(SEED)

    refs_dir = os.path.join(OUTPUT_DIR, "references")
    fastq_dir = os.path.join(OUTPUT_DIR, "fastq")
    os.makedirs(refs_dir, exist_ok=True)
    os.makedirs(fastq_dir, exist_ok=True)

    # ── Reference sequences + SNP definitions (shared across samples) ──
    locus_data = {}
    for idx, (locus_name, length) in enumerate(LOCI.items()):
        ref_seq = "".join(rng.choice(BASES) for _ in range(length))
        positions = snp_positions(length, idx)
        alt_alleles = {pos: TRANSITIONS[ref_seq[pos]] for pos in positions}
        locus_data[locus_name] = {
            "ref_seq": ref_seq,
            "snp_positions": positions,
            "alt_alleles": alt_alleles,
        }

        # Per-locus FASTA: one record per sample (all share the reference seq)
        fasta_path = os.path.join(refs_dir, f"{locus_name}.fasta")
        with open(fasta_path, "w") as f:
            for sample in SAMPLES:
                f.write(f">{sample['id']}\n{ref_seq}\n")
        print(f"  Wrote {fasta_path}")

    # ── Per-sample paired-end FASTQ ────────────────────────────────
    for sample in SAMPLES:
        sid, ploidy, alt_dosage = sample["id"], sample["ploidy"], sample["alt_dosage"]
        alt_haplotypes = set(range(ploidy - alt_dosage, ploidy))  # which haps carry alt

        r1_path = os.path.join(fastq_dir, f"{sid}_R1.fastq.gz")
        r2_path = os.path.join(fastq_dir, f"{sid}_R2.fastq.gz")
        read_num = 0
        with gzip.open(r1_path, "wt") as r1_fh, gzip.open(r2_path, "wt") as r2_fh:
            for locus_name, data in locus_data.items():
                ref_seq = data["ref_seq"]
                length = len(ref_seq)
                max_start = max(0, length - FRAGMENT_SIZE)

                for pair_idx in range(NUM_READ_PAIRS_PER_LOCUS):
                    read_num += 1
                    frag_start = rng.randint(0, max_start)
                    fragment = list(ref_seq[frag_start:frag_start + FRAGMENT_SIZE])

                    haplotype = pair_idx % ploidy
                    if haplotype in alt_haplotypes:
                        for snp_pos in data["snp_positions"]:
                            local = snp_pos - frag_start
                            if 0 <= local < len(fragment):
                                fragment[local] = data["alt_alleles"][snp_pos]
                    fragment = "".join(fragment)

                    r1_seq = fragment[:READ_LENGTH]
                    r2_seq = reverse_complement(fragment[-READ_LENGTH:])
                    qual = QUAL_CHAR * READ_LENGTH
                    name = f"@{sid}:{locus_name}:{read_num}"
                    r1_fh.write(f"{name}/1\n{r1_seq}\n+\n{qual}\n")
                    r2_fh.write(f"{name}/2\n{r2_seq}\n+\n{qual}\n")

        print(f"  Wrote {r1_path} / {r2_path} ({read_num} read pairs)")

    # ── Samplesheet (absolute FASTQ paths) ─────────────────────────
    samplesheet_path = os.path.join(OUTPUT_DIR, "samplesheet.csv")
    with open(samplesheet_path, "w") as f:
        f.write("sample,fastq_1,fastq_2,species,ploidy\n")
        for sample in SAMPLES:
            sid = sample["id"]
            r1 = os.path.join(fastq_dir, f"{sid}_R1.fastq.gz")
            r2 = os.path.join(fastq_dir, f"{sid}_R2.fastq.gz")
            f.write(f"{sid},{r1},{r2},{SPECIES},{sample['ploidy']}\n")
    print(f"  Wrote {samplesheet_path}")

    # ── Summary ────────────────────────────────────────────────────
    print("\nTest data summary:")
    summary = ", ".join(
        f"{s['id']}(p{s['ploidy']},alt{s['alt_dosage']})" for s in SAMPLES
    )
    print(f"  Samples: {summary}")
    print(f"  Loci: {len(LOCI)} ({', '.join(f'{k}:{v}bp' for k, v in LOCI.items())})")
    for locus_name, data in locus_data.items():
        snps = ", ".join(
            f"{pos}({data['ref_seq'][pos]}->{data['alt_alleles'][pos]})"
            for pos in data["snp_positions"]
        )
        print(f"  {locus_name} SNPs: {snps}")


if __name__ == "__main__":
    main()
