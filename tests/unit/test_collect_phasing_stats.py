"""Unit tests for bin/collect_phasing_stats.py.

These lock down the round-8 rewrite: the per-(sample, locus) O(N^2) rescans were
replaced with single-pass precomputes, so correctness of those precomputes is the
key thing to verify.
"""
import gzip

import collect_phasing_stats as cps


def test_build_allele_counts_single_pass(tmp_path):
    (tmp_path / "sample1.locus1.phased.fasta").write_text(
        ">sample1__locus1__1\nACGT\n>sample1__locus1__2\nACGA\n"
    )
    (tmp_path / "sample1.locus2.phased.fasta").write_text(
        ">sample1__locus2__1\nTTTT\n"
    )
    (tmp_path / "sample2.locus1.phased.fasta").write_text(
        ">sample2__locus1__1\nGGGG\n>sample2__locus1__2\nGGGA\n"
        ">sample2__locus1__3\nGGAA\n>sample2__locus1__4\nGAAA\n"
    )
    counts = cps.build_allele_counts(str(tmp_path))
    assert counts[("sample1", "locus1")] == 2
    assert counts[("sample1", "locus2")] == 1
    assert counts[("sample2", "locus1")] == 4
    assert ("sample2", "locus2") not in counts


def test_parse_ref_lengths(tmp_path):
    ref = tmp_path / "sample1.ref.fasta"
    ref.write_text(">locus1\nACGTAC\n>locus2\nTT\nTT\n")   # 6 and 4 (multi-line)
    lengths = cps.parse_ref_lengths(str(ref))
    assert lengths == {"locus1": 6, "locus2": 4}


def test_parse_ref_lengths_missing_file_returns_empty(tmp_path):
    assert cps.parse_ref_lengths(str(tmp_path / "nope.fasta")) == {}


def test_count_variants_by_locus_plain_and_gz(tmp_path):
    body = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"
        "locusA\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "locusA\t20\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n"
        "locusA\t30\t.\tG\tA\t50\tlowqual\t.\tGT\t0/1\n"   # not PASS
        "locusB\t10\t.\tT\tC\t50\tPASS\t.\tGT\t0/1\n"
    )
    plain = tmp_path / "s.vcf"
    plain.write_text(body)
    counts = cps.count_variants_by_locus(str(plain))
    assert counts["locusA"] == 2      # PASS only
    assert counts["locusB"] == 1

    gz = tmp_path / "s.vcf.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(body)
    assert dict(cps.count_variants_by_locus(str(gz))) == {"locusA": 2, "locusB": 1}


def test_count_variants_by_locus_missing_file(tmp_path):
    assert cps.count_variants_by_locus(str(tmp_path / "nope.vcf")) == {}


def test_parse_phase_blocks(tmp_path):
    p = tmp_path / "b.phase.out"
    # one block of length 3 (declared in the header via "len: 3")
    p.write_text("#BLOCK: offset: 1 len: 3\n1\t0\t1\n2\t1\t0\n3\t0\t1\n********\n")
    n_blocks, longest = cps.parse_phase_blocks(str(p))
    assert n_blocks == 1
    assert longest == 3


def test_parse_phase_blocks_empty(tmp_path):
    empty = tmp_path / "e.phase.out"
    empty.write_text("")
    assert cps.parse_phase_blocks(str(empty)) == (0, 0)


def test_find_sample_vcf_prefers_exact_name(tmp_path):
    (tmp_path / "sample1.snps.biallelic.filtered.vcf.gz").write_text("x")
    got = cps.find_sample_vcf(str(tmp_path), "sample1")
    assert got.endswith("sample1.snps.biallelic.filtered.vcf.gz")
