#!/usr/bin/env python3
"""
Collect phasing statistics across all samples and loci.

Generates per-individual summary files and an aggregate
averagePhasingStats.txt file.

Replaces PATE.pl lines 1100-1391 (summary statistics logic).
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict


def count_pass_variants(vcf_path):
    """Count PASS variants in a VCF file."""
    count = 0
    if not os.path.exists(vcf_path):
        return 0
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 7 and fields[6] == 'PASS':
                count += 1
    return count


def parse_phase_blocks(phase_path):
    """Parse H-PoPG .phase.out file for block statistics."""
    n_blocks = 0
    longest_block = 0
    current_block_len = 0

    if not os.path.exists(phase_path):
        return 0, 0

    with open(phase_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#BLOCK:'):
                if current_block_len > longest_block:
                    longest_block = current_block_len
                n_blocks += 1
                current_block_len = 0
                # Also try to parse block length from header
                match = re.search(r'len:\s+(\d+)', line)
                if match:
                    current_block_len = int(match.group(1))
            elif line.startswith('*'):
                if current_block_len > longest_block:
                    longest_block = current_block_len
            elif line and line[0].isdigit():
                # Count lines in block if not parsed from header
                pass

    if current_block_len > longest_block:
        longest_block = current_block_len

    return n_blocks, longest_block


def get_ref_length(ref_fasta, locus):
    """Get reference sequence length for a locus."""
    current_header = None
    current_len = 0

    if not os.path.exists(ref_fasta):
        return 0

    with open(ref_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header == locus:
                    return current_len
                current_header = line[1:].split()[0]
                current_len = 0
            elif line and current_header == locus:
                current_len += len(line)

    if current_header == locus:
        return current_len
    return 0


def count_alleles_in_phased(phased_dir, sample_id, locus):
    """Count the number of alleles for a sample in a phased FASTA."""
    count = 0
    pattern = f"{sample_id}__{locus}__"

    for f in os.listdir(phased_dir):
        if f.endswith('.phased.fasta'):
            filepath = os.path.join(phased_dir, f)
            with open(filepath, 'r') as fh:
                for line in fh:
                    if line.startswith('>') and pattern in line:
                        count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description='Collect phasing statistics')
    parser.add_argument('--phased_dir', required=True, help='Directory with phased FASTA files')
    parser.add_argument('--vcf_dir', required=True, help='Directory with per-locus VCF files')
    parser.add_argument('--phase_dir', required=True, help='Directory with H-PoPG .phase.out files')
    parser.add_argument('--reference_dir', required=True, help='Directory with per-sample reference FASTAs')
    parser.add_argument('--ploidy_json', required=True, help='JSON list of {id, ploidy} maps')
    parser.add_argument('--output_dir', required=True, help='Output directory for stats files')
    args = parser.parse_args()

    # Parse ploidy information
    ploidy_data = json.loads(args.ploidy_json)
    sample_ploidies = {}
    for entry in ploidy_data:
        sample_ploidies[entry['id']] = int(entry['ploidy'])

    max_ploidy = max(sample_ploidies.values()) if sample_ploidies else 2

    # Discover loci from phased FASTA files
    loci = set()
    for f in os.listdir(args.phased_dir):
        if f.endswith('.phased.fasta'):
            with open(os.path.join(args.phased_dir, f), 'r') as fh:
                for line in fh:
                    if line.startswith('>'):
                        parts = line[1:].strip().split('__')
                        if len(parts) >= 2:
                            loci.add(parts[1])

    # Also discover from VCF files
    for f in os.listdir(args.vcf_dir):
        if f.endswith('.vcf'):
            # Extract locus from filename pattern: sample.snps.biallelic.LOCUS.vcf
            parts = f.split('.')
            for i, p in enumerate(parts):
                if p == 'biallelic' and i + 1 < len(parts) - 1:
                    loci.add(parts[i + 1])

    loci = sorted(loci)

    # Collect stats per sample
    global_stats = {}
    global_counts = {}

    for sample_id, ploidy in sorted(sample_ploidies.items()):
        stats_file = os.path.join(args.output_dir, f"{sample_id}.phasingSummary.txt")

        # Initialize global stats
        global_stats[sample_id] = {
            'nloci': 0, 'ninvarloci': 0, 'nvarloci': 0, 'nphaseloci': 0,
            'loclen': 0, 'nvar': 0, 'het': 0, 'nblocks': 0, 'lbl': 0
        }
        global_counts[sample_id] = defaultdict(int)

        # Find reference FASTA for this sample
        ref_fasta = None
        for f in os.listdir(args.reference_dir):
            if f.startswith(sample_id) and f.endswith('.ref.fasta'):
                ref_fasta = os.path.join(args.reference_dir, f)
                break

        with open(stats_file, 'w') as out:
            # Header
            header_parts = ['LOCUS', 'LENGTH', 'NVAR', 'HET', 'NBLOCKS', 'LONGESTBL']
            for i in range(1, max_ploidy + 1):
                header_parts.append(str(i))
            out.write('\t'.join(header_parts) + '\n')

            for locus in loci:
                ref_len = get_ref_length(ref_fasta, locus) if ref_fasta else 0

                # Count variants from VCF
                vcf_pattern = f"{sample_id}.snps.biallelic.filtered.{locus}.vcf"
                vcf_path = os.path.join(args.vcf_dir, vcf_pattern)
                if not os.path.exists(vcf_path):
                    # Try alternative naming
                    for f in os.listdir(args.vcf_dir):
                        if sample_id in f and locus in f and f.endswith('.vcf'):
                            vcf_path = os.path.join(args.vcf_dir, f)
                            break

                nvar = count_pass_variants(vcf_path)

                # Parse phase blocks
                phase_pattern = f"{sample_id}.{locus}.phase.out"
                phase_path = os.path.join(args.phase_dir, phase_pattern)
                nblocks, lbl = parse_phase_blocks(phase_path)

                # Calculate heterozygosity
                het = nvar / ref_len if ref_len > 0 else 0

                # Count alleles
                n_alleles = count_alleles_in_phased(args.phased_dir, sample_id, locus)

                # Write per-locus stats
                row = [locus, str(ref_len), str(nvar), f"{het:.6f}", str(nblocks), str(lbl)]
                for i in range(1, max_ploidy + 1):
                    if n_alleles == i:
                        row.append('1')
                        global_counts[sample_id][i] += 1
                    elif n_alleles > 0:
                        row.append('0')
                    else:
                        row.append('NA')
                out.write('\t'.join(row) + '\n')

                # Update global stats
                if n_alleles > 0:
                    global_stats[sample_id]['nloci'] += 1
                    global_stats[sample_id]['loclen'] += ref_len

                    if nvar > 0 and nblocks > 0:
                        global_stats[sample_id]['nvar'] += nvar
                        global_stats[sample_id]['het'] += het
                        global_stats[sample_id]['nblocks'] += nblocks
                        global_stats[sample_id]['lbl'] += lbl
                        global_stats[sample_id]['nphaseloci'] += 1
                        global_stats[sample_id]['nvarloci'] += 1
                    elif nvar > 0 and nblocks == 0:
                        global_stats[sample_id]['nvar'] += nvar
                        global_stats[sample_id]['het'] += het
                        global_stats[sample_id]['nvarloci'] += 1
                    elif nvar == 0:
                        global_stats[sample_id]['ninvarloci'] += 1

    # Write aggregate stats
    avg_file = os.path.join(args.output_dir, 'averagePhasingStats.txt')
    with open(avg_file, 'w') as out:
        header = ['INDIVIDUAL', 'NLOCI', 'NVARLOCI', 'NINVLOCI', 'NPHASELOCI',
                  'AVG_LENGTH', 'AVG_NVAR', 'AVG_HET', 'AVG_NBLOCKS', 'AVG_LONGESTBL']
        for i in range(1, max_ploidy + 1):
            header.append(str(i))
        out.write('\t'.join(header) + '\n')

        for sample_id in sorted(sample_ploidies.keys()):
            gs = global_stats[sample_id]
            ploidy = sample_ploidies[sample_id]

            if gs['nloci'] > 0:
                avg_len = gs['loclen'] / gs['nloci']
                avg_nvar = gs['nvar'] / gs['nloci']
                avg_het = gs['het'] / gs['nloci']
                avg_nblocks = gs['nblocks'] / gs['nphaseloci'] if gs['nphaseloci'] > 0 else 0
                avg_lbl = gs['lbl'] / gs['nphaseloci'] if gs['nphaseloci'] > 0 else 0

                row = [sample_id, str(gs['nloci']), str(gs['nvarloci']),
                       str(gs['ninvarloci']), str(gs['nphaseloci']),
                       f"{avg_len:.2f}", f"{avg_nvar:.4f}", f"{avg_het:.6f}",
                       f"{avg_nblocks:.2f}", f"{avg_lbl:.2f}"]
            else:
                row = [sample_id, '0', '0', '0', '0', '0', '0', '0', '0', '0']

            for i in range(1, max_ploidy + 1):
                if i <= ploidy:
                    row.append(str(global_counts[sample_id].get(i, 0)))
                else:
                    row.append('NA')

            out.write('\t'.join(row) + '\n')

    print(f"Statistics collected for {len(sample_ploidies)} samples across {len(loci)} loci")


if __name__ == '__main__':
    main()
