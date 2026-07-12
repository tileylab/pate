#!/usr/bin/env python3
"""
Build phased consensus sequences from H-PoPG output and VCF.

For a given sample and locus, parses the H-PoPG phase output file
to determine phased haplotype blocks, reads the VCF for allele
positions, and applies phased alleles to the reference (or IUPAC)
sequence to produce:
  - PHASED FASTA: phased haplotype sequences
  - GENOTYPE FASTA: unphased IUPAC-coded sequence
  - PICKONE FASTA: one randomly chosen allele

Replaces PATE.pl lines 581-1060 (alleles mode logic).
"""

import argparse
import gzip
import os
import random
import sys


def parse_fasta_single(filepath):
    """Parse a FASTA file, return dict of {locus: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                # Handle GATK FastaAlternateReferenceMaker header format: >1 locusname:1-100
                parts = line[1:].split()
                if len(parts) >= 2 and ':' in parts[1]:
                    current_header = parts[1].split(':')[0]
                else:
                    current_header = parts[0]
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_header is not None:
        sequences[current_header] = ''.join(current_seq)

    return sequences


def parse_vcf(vcf_path, locus=None):
    """Parse VCF and return PASS variants (pos, ref, alt) for one locus.

    The per-sample VCF contains variants for every locus, and each locus is a
    separate contig whose POS restarts at 1. Restrict to ``locus`` (CHROM) so
    positions map onto the correct reference and stay aligned with the per-locus
    order H-PoPG phased against; otherwise variants from other loci collide onto
    this locus's coordinates and mask most of it to N.
    """
    variants = []
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, filt = fields[:7]
            if filt == 'PASS' and (locus is None or chrom == locus):
                variants.append({
                    'pos': int(pos),
                    'ref': ref,
                    'alt': alt
                })
    return variants


def parse_phase_output(phase_path):
    """
    Parse H-PoPG .phase.out file.

    Returns:
        n_chunks: number of phasing blocks
        chunk_lengths: dict of {chunk_id: length}
        max_chunk: chunk_id of the longest block
        snp_matrix: dict of {snp_index: {allele_index: value}}
        chunks: dict of {chunk_id: [snp_indices]}
        ploidy: detected ploidy from the phase file
    """
    chunks = {}
    chunk_lengths = {}
    snp_matrix = {}
    n_chunks = 0
    max_chunk = 0
    max_chunk_len = 0
    ploidy = 0

    if not os.path.exists(phase_path):
        return 0, {}, 0, {}, {}, 0

    with open(phase_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#BLOCK:'):
                n_chunks += 1
                chunk_lengths[n_chunks] = 0
            elif line.startswith('*'):
                if chunk_lengths.get(n_chunks, 0) > max_chunk_len:
                    max_chunk_len = chunk_lengths[n_chunks]
                    max_chunk = n_chunks
            elif line and line[0].isdigit():
                parts = line.split()
                snp_idx = int(parts[0])
                alleles = parts[1:]
                chunk_lengths[n_chunks] = chunk_lengths.get(n_chunks, 0) + 1

                if ploidy == 0:
                    ploidy = len(alleles)

                if n_chunks not in chunks:
                    chunks[n_chunks] = []
                chunks[n_chunks].append(snp_idx)

                snp_matrix[snp_idx] = {}
                for i, allele in enumerate(alleles):
                    snp_matrix[snp_idx][i + 1] = allele

    # Check last chunk
    if n_chunks > 0 and chunk_lengths.get(n_chunks, 0) > max_chunk_len:
        max_chunk = n_chunks

    return n_chunks, chunk_lengths, max_chunk, snp_matrix, chunks, ploidy


def build_phased_sequences(ref_seq, variants, snp_matrix, chunks, max_chunk,
                           ploidy, genotype_mode, iupac_seq=None):
    """Build phased haplotype sequences."""
    # Map variant positions to indices (1-based)
    pos_to_idx = {}
    idx_to_pos = {}
    alleles = {}
    for i, var in enumerate(variants):
        idx = i + 1
        pos_to_idx[var['pos']] = idx
        idx_to_pos[idx] = var['pos']
        alleles[var['pos']] = {0: var['ref'], 1: var['alt']}

    # Get set of SNP indices in the longest block
    max_chunk_snps = set(chunks.get(max_chunk, []))

    phased_seqs = []
    for allele_num in range(1, ploidy + 1):
        base_seq = ref_seq if genotype_mode == 'consensus' else (iupac_seq or ref_seq)
        seq_chars = list(base_seq)
        phased_chars = []

        for j, char in enumerate(seq_chars):
            pos = j + 1
            if pos not in pos_to_idx:
                phased_chars.append(char)
            else:
                snp_idx = pos_to_idx[pos]
                if snp_idx in max_chunk_snps:
                    # This SNP is in the longest phased block
                    if snp_idx in snp_matrix and allele_num in snp_matrix[snp_idx]:
                        val = snp_matrix[snp_idx][allele_num]
                        if val.isdigit():
                            val_int = int(val)
                            if val_int in alleles[pos]:
                                phased_chars.append(alleles[pos][val_int])
                            else:
                                phased_chars.append('N')
                        elif val == '-':
                            phased_chars.append('N')
                        else:
                            phased_chars.append('N')
                    else:
                        phased_chars.append('N')
                else:
                    # SNP not in longest block - use N for consensus, keep IUPAC for iupac mode
                    if genotype_mode == 'consensus':
                        phased_chars.append('N')
                    else:
                        phased_chars.append(char)

        header = f"{allele_num}"
        phased_seqs.append((header, ''.join(phased_chars)))

    return phased_seqs


def main():
    parser = argparse.ArgumentParser(description='Build phased consensus sequences')
    parser.add_argument('--phase_file', required=True, help='H-PoPG .phase.out file')
    parser.add_argument('--vcf', required=True, help='Per-locus filtered VCF file')
    parser.add_argument('--reference_fasta', required=True, help='Per-sample reference FASTA')
    parser.add_argument('--iupac_fasta', required=True, help='IUPAC consensus FASTA')
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--locus', required=True, help='Locus name')
    parser.add_argument('--ploidy', required=True, type=int, help='Sample ploidy')
    parser.add_argument('--genotype_mode', default='consensus', choices=['consensus', 'iupac'],
                        help='Genotype mode: consensus or iupac')
    parser.add_argument('--unique_only', action='store_true', help='Output only unique haplotypes')
    parser.add_argument('--output_expected_dosage', action='store_true',
                        help='Output expected dosage copies for invariant loci')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    args = parser.parse_args()

    # Parse reference sequences
    ref_seqs = parse_fasta_single(args.reference_fasta)
    iupac_seqs = parse_fasta_single(args.iupac_fasta)

    # Get reference sequence for this locus
    ref_seq = ref_seqs.get(args.locus, '')
    iupac_seq = iupac_seqs.get(args.locus, '')

    if not ref_seq:
        print(f"WARNING: Locus {args.locus} not found in reference for sample {args.sample_id}",
              file=sys.stderr)
        # Write empty output files
        for suffix in ['phased', 'genotype', 'pickone']:
            with open(f"{args.output_prefix}.{suffix}.fasta", 'w') as f:
                pass
        return

    # Parse VCF — only this locus's variants (POS is 1-based within the locus contig)
    variants = parse_vcf(args.vcf, args.locus)

    # Parse phase output
    n_chunks, chunk_lengths, max_chunk, snp_matrix, chunks, detected_ploidy = \
        parse_phase_output(args.phase_file)

    # Build phased sequences
    phased_fasta = f"{args.output_prefix}.phased.fasta"
    genotype_fasta = f"{args.output_prefix}.genotype.fasta"
    pickone_fasta = f"{args.output_prefix}.pickone.fasta"

    with open(phased_fasta, 'w') as pf, \
         open(pickone_fasta, 'w') as pkf:

        if len(variants) > 0 and n_chunks > 0:
            # There are phased variants
            phased_seqs = build_phased_sequences(
                ref_seq, variants, snp_matrix, chunks, max_chunk,
                args.ploidy, args.genotype_mode, iupac_seq
            )

            # Apply unique_only filter
            if args.unique_only:
                seen_seqs = set()
                unique_seqs = []
                for header, seq in phased_seqs:
                    if seq not in seen_seqs:
                        seen_seqs.add(seq)
                        unique_seqs.append((header, seq))
                phased_seqs = unique_seqs

            # Write phased sequences
            for header, seq in phased_seqs:
                full_header = f"{args.sample_id}__{args.locus}__{header}"
                pf.write(f">{full_header}\n{seq}\n")

            # Write pickone (random allele)
            if phased_seqs:
                rand_idx = random.randint(0, len(phased_seqs) - 1)
                header, seq = phased_seqs[rand_idx]
                full_header = f"{args.sample_id}__{args.locus}__{header}"
                pkf.write(f">{full_header}\n{seq}\n")

        else:
            # No variants or no phasing - output reference/IUPAC sequence
            base_seq = ref_seq if args.genotype_mode == 'consensus' else (iupac_seq or ref_seq)

            if args.output_expected_dosage:
                for i in range(1, args.ploidy + 1):
                    header = f"{args.sample_id}__{args.locus}__{i}"
                    pf.write(f">{header}\n{base_seq}\n")
            else:
                header = f"{args.sample_id}__{args.locus}__REF"
                pf.write(f">{header}\n{base_seq}\n")

            ref_header = f"{args.sample_id}__{args.locus}__REF"
            pkf.write(f">{ref_header}\n{base_seq}\n")

    # Write genotype FASTA (unphased IUPAC-coded sequence)
    with open(genotype_fasta, 'w') as gf:
        genotype_seq = iupac_seq if iupac_seq else ref_seq
        gf.write(f">{args.sample_id}\n{genotype_seq}\n")


if __name__ == '__main__':
    main()
