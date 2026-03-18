#!/usr/bin/env python3
"""
Prepare a per-sample concatenated reference FASTA from per-locus FASTA files.

For each locus FASTA in the references directory, extracts the sequence
matching the given sample ID (or species name) and writes all matching
sequences into a single concatenated reference FASTA file.

Replaces PATE.pl lines 278-333 (reference preparation logic).
"""

import argparse
import os
import sys


def parse_fasta(filepath):
    """Parse a FASTA file and return dict of {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                # Strip HybPiper-style suffixes like -L001
                if '-L' in current_header:
                    current_header = current_header.rsplit('-L', 1)[0]
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_header is not None:
        sequences[current_header] = ''.join(current_seq)

    return sequences


def main():
    parser = argparse.ArgumentParser(description='Prepare per-sample reference from per-locus FASTAs')
    parser.add_argument('--references_dir', required=True, help='Directory containing per-locus FASTA files')
    parser.add_argument('--sample_id', required=True, help='Sample ID to extract from each locus FASTA')
    parser.add_argument('--species', required=False, default=None, help='Species name (alternative match key)')
    parser.add_argument('--output', required=True, help='Output concatenated reference FASTA path')
    args = parser.parse_args()

    refs_dir = args.references_dir
    sample_id = args.sample_id
    species = args.species
    output = args.output

    locus_count = 0

    # Find all FASTA files in the references directory
    fasta_files = sorted([
        f for f in os.listdir(refs_dir)
        if f.endswith('.fasta') or f.endswith('.fa') or f.endswith('.fna')
    ])

    if not fasta_files:
        print(f"ERROR: No FASTA files found in {refs_dir}", file=sys.stderr)
        sys.exit(1)

    with open(output, 'w') as out_fh:
        for fasta_file in fasta_files:
            locus = os.path.splitext(fasta_file)[0]
            filepath = os.path.join(refs_dir, fasta_file)
            sequences = parse_fasta(filepath)

            # Try matching by sample_id first, then by species
            seq = None
            if sample_id in sequences:
                seq = sequences[sample_id]
            elif species and species in sequences:
                seq = sequences[species]

            if seq is not None:
                out_fh.write(f">{locus}\n{seq}\n")
                locus_count += 1

    if locus_count == 0:
        print(f"WARNING: No matching sequences found for sample '{sample_id}' "
              f"(species: '{species}') in {refs_dir}", file=sys.stderr)
    else:
        print(f"Prepared reference with {locus_count} loci for sample {sample_id}")


if __name__ == '__main__':
    main()
