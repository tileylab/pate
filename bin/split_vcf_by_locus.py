#!/usr/bin/env python3
"""
Extract per-locus VCF regions from a compressed, indexed VCF.

This is a simple wrapper that uses tabix to extract variants
for a specific genomic region (locus/contig) from a bgzipped VCF.

Replaces the inline tabix commands in PATE.pl.
"""

import argparse
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description='Extract per-locus VCF from bgzipped VCF')
    parser.add_argument('--vcf', required=True, help='Input bgzipped VCF file (.vcf.gz)')
    parser.add_argument('--locus', required=True, help='Locus/contig name to extract')
    parser.add_argument('--output', required=True, help='Output VCF file path')
    args = parser.parse_args()

    # Extract region with tabix (include header with -h flag)
    try:
        result = subprocess.run(
            ['tabix', '-h', args.vcf, args.locus],
            capture_output=True, text=True, check=True
        )
        with open(args.output, 'w') as out:
            out.write(result.stdout)
    except subprocess.CalledProcessError as e:
        # If tabix fails (e.g. no variants in region), write just the header
        result = subprocess.run(
            ['tabix', '-H', args.vcf],
            capture_output=True, text=True
        )
        with open(args.output, 'w') as out:
            out.write(result.stdout)
    except FileNotFoundError:
        print("ERROR: tabix not found in PATH", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
