#!/usr/bin/env python3
"""
Collect phasing statistics across all samples and loci.

Generates per-individual summary files and an aggregate
averagePhasingStats.txt file.

Replaces PATE.pl lines 1100-1391 (summary statistics logic).
"""

import argparse
import gzip
import json
import os
import re
import sys
from collections import defaultdict


def find_sample_vcf(vcf_dir, sample_id):
    """Locate a sample's filtered VCF (gzipped or not) in the collected dir."""
    for name in (f"{sample_id}.snps.biallelic.filtered.vcf.gz",
                 f"{sample_id}.snps.biallelic.filtered.vcf"):
        path = os.path.join(vcf_dir, name)
        if os.path.exists(path):
            return path
    # Fallback: any per-sample filtered VCF for this sample
    for f in sorted(os.listdir(vcf_dir)):
        if f.startswith(sample_id + '.') and 'filtered' in f and (f.endswith('.vcf') or f.endswith('.vcf.gz')):
            return os.path.join(vcf_dir, f)
    return None


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


def build_allele_counts(phased_dir):
    """Scan the phased dir once, returning {(sample_id, locus): allele_count}.

    Headers are `>sample__locus__hap`. Doing this in a single pass avoids an
    O(files^2) rescan of every phased FASTA for every (sample, locus).
    """
    counts = defaultdict(int)
    for f in os.listdir(phased_dir):
        if f.endswith('.phased.fasta'):
            with open(os.path.join(phased_dir, f), 'r') as fh:
                for line in fh:
                    if line.startswith('>'):
                        parts = line[1:].strip().split('__')
                        if len(parts) >= 2:
                            counts[(parts[0], parts[1])] += 1
    return counts


def parse_ref_lengths(ref_fasta):
    """Parse a per-sample reference FASTA once, returning {locus: length}."""
    lengths = {}
    if not ref_fasta or not os.path.exists(ref_fasta):
        return lengths
    header = None
    length = 0
    with open(ref_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    lengths[header] = length
                header = line[1:].split()[0]
                length = 0
            elif line:
                length += len(line)
    if header is not None:
        lengths[header] = length
    return lengths


def count_variants_by_locus(vcf_path):
    """Parse a per-sample VCF once, returning {locus(CHROM): PASS variant count}."""
    counts = defaultdict(int)
    if not vcf_path or not os.path.exists(vcf_path):
        return counts
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 7 and fields[6] == 'PASS':
                counts[fields[0]] += 1
    return counts


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

    # Discover loci from phased FASTA headers (sample__locus__hap)
    loci = set()
    for f in os.listdir(args.phased_dir):
        if f.endswith('.phased.fasta'):
            with open(os.path.join(args.phased_dir, f), 'r') as fh:
                for line in fh:
                    if line.startswith('>'):
                        parts = line[1:].strip().split('__')
                        if len(parts) >= 2:
                            loci.add(parts[1])

    loci = sorted(loci)

    # Precompute allele counts once for all (sample, locus) — avoids O(files^2).
    allele_counts = build_allele_counts(args.phased_dir)

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

        # Per-sample lookups, each parsed once (not per locus)
        ref_lengths = parse_ref_lengths(ref_fasta)
        nvar_by_locus = count_variants_by_locus(find_sample_vcf(args.vcf_dir, sample_id))

        with open(stats_file, 'w') as out:
            # Header
            header_parts = ['LOCUS', 'LENGTH', 'NVAR', 'HET', 'NBLOCKS', 'LONGESTBL']
            for i in range(1, max_ploidy + 1):
                header_parts.append(str(i))
            out.write('\t'.join(header_parts) + '\n')

            for locus in loci:
                ref_len = ref_lengths.get(locus, 0)

                # This locus's PASS variants from the per-sample VCF (by CHROM)
                nvar = nvar_by_locus.get(locus, 0)

                # Parse phase blocks
                phase_pattern = f"{sample_id}.{locus}.phase.out"
                phase_path = os.path.join(args.phase_dir, phase_pattern)
                nblocks, lbl = parse_phase_blocks(phase_path)

                # Calculate heterozygosity
                het = nvar / ref_len if ref_len > 0 else 0

                # Allele count (precomputed)
                n_alleles = allele_counts.get((sample_id, locus), 0)

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
