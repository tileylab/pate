"""Shared pytest fixtures for the bin/ script unit tests.

The pipeline's helper scripts live in ``bin/`` as standalone executables (their
``main()`` is guarded by ``if __name__ == '__main__'``), so they import cleanly
as modules once ``bin/`` is on ``sys.path``.
"""
import os
import sys

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN_DIR = os.path.join(REPO_ROOT, "bin")

if BIN_DIR not in sys.path:
    sys.path.insert(0, BIN_DIR)


@pytest.fixture
def multi_locus_vcf(tmp_path):
    """A tiny VCF with PASS variants on two loci at overlapping positions.

    Each locus is its own contig with POS restarting at 1, so positions collide
    across loci — exactly the situation that must be filtered by CHROM.
    """
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
        # locusA
        "locusA\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
        "locusA\t20\t.\tC\tT\t50\tPASS\t.\tGT\t0/1",
        "locusA\t30\t.\tG\tA\t50\tmyfilter\t.\tGT\t0/1",   # not PASS
        # locusB — same POS values as locusA
        "locusB\t10\t.\tT\tC\t50\tPASS\t.\tGT\t0/1",
        "locusB\t20\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
        "locusB\t25\t.\tC\tA\t50\tPASS\t.\tGT\t0/1",
    ]
    p = tmp_path / "sample1.snps.biallelic.filtered.vcf"
    p.write_text("\n".join(lines) + "\n")
    return str(p)
