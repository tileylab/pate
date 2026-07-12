"""Unit tests for bin/build_phased_consensus.py.

The headline test is the round-8 regression guard: ``parse_vcf`` must restrict
variants to the requested locus, otherwise variants from other loci (whose POS
restarts at 1 per contig) collide onto the current locus and mask most of it to N.
"""
import gzip

import build_phased_consensus as bpc


# ── parse_vcf: the cross-locus regression guard ─────────────────────────────
def test_parse_vcf_filters_to_locus(multi_locus_vcf):
    a = bpc.parse_vcf(multi_locus_vcf, "locusA")
    assert [v["pos"] for v in a] == [10, 20]          # only locusA PASS records
    assert all(v["ref"] and v["alt"] for v in a)

    b = bpc.parse_vcf(multi_locus_vcf, "locusB")
    assert [v["pos"] for v in b] == [10, 20, 25]       # locusB PASS records


def test_parse_vcf_without_locus_returns_all_pass(multi_locus_vcf):
    allv = bpc.parse_vcf(multi_locus_vcf)
    # 2 PASS on locusA + 3 PASS on locusB = 5 (the non-PASS record is excluded)
    assert len(allv) == 5


def test_parse_vcf_reads_gzip(tmp_path):
    gz = tmp_path / "s.vcf.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n")
        fh.write("locusA\t5\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n")
    assert bpc.parse_vcf(str(gz), "locusA") == [{"pos": 5, "ref": "A", "alt": "T"}]


# ── build_phased_sequences: correctness + not-mostly-N guard ─────────────────
def _phase_two_hap():
    """ref + two het variants fully phased into one block, ploidy 2."""
    ref = "ACGTACGTAC"                       # 10 bp
    variants = [
        {"pos": 3, "ref": "G", "alt": "C"},  # ref[2] == 'G'
        {"pos": 7, "ref": "G", "alt": "T"},  # ref[6] == 'G'
    ]
    # phase.out: snp1 alleles "0 1", snp2 alleles "1 0" (allele_num is 1-based)
    snp_matrix = {1: {1: "0", 2: "1"}, 2: {1: "1", 2: "0"}}
    chunks = {1: [1, 2]}
    return ref, variants, snp_matrix, chunks


def test_build_phased_sequences_copies_reference_and_applies_alleles():
    ref, variants, snp_matrix, chunks = _phase_two_hap()
    seqs = bpc.build_phased_sequences(
        ref, variants, snp_matrix, chunks, max_chunk=1,
        ploidy=2, genotype_mode="consensus",
    )
    assert [h for h, _ in seqs] == ["1", "2"]
    hap1 = seqs[0][1]
    hap2 = seqs[1][1]

    # length preserved, and (this locus) fully phased -> no N masking
    assert len(hap1) == len(ref) == 10
    assert "N" not in hap1 and "N" not in hap2

    # non-variant positions equal the reference in both haplotypes
    variant_idx = {v["pos"] - 1 for v in variants}   # 0-based
    for i in range(len(ref)):
        if i not in variant_idx:
            assert hap1[i] == ref[i]
            assert hap2[i] == ref[i]

    # variant positions carry the phased allele: snp1 -> hap1=ref/hap2=alt, snp2 -> hap1=alt/hap2=ref
    assert hap1[2] == "G" and hap2[2] == "C"   # pos 3
    assert hap1[6] == "T" and hap2[6] == "G"   # pos 7


def test_build_phased_sequences_not_mostly_n_for_normal_locus():
    """Regression guard for the round-8 symptom (haplotypes were ~90% N)."""
    ref, variants, snp_matrix, chunks = _phase_two_hap()
    seqs = bpc.build_phased_sequences(
        ref, variants, snp_matrix, chunks, max_chunk=1,
        ploidy=2, genotype_mode="consensus",
    )
    for _, seq in seqs:
        assert seq.count("N") / len(seq) < 0.5


# ── parse_fasta_single: GATK vs plain headers ───────────────────────────────
def test_parse_fasta_single_handles_gatk_and_plain_headers(tmp_path):
    p = tmp_path / "iupac.fasta"
    p.write_text(
        ">1 4471:1-8\nACGTACGT\n"     # GATK FastaAlternateReferenceMaker style
        ">locus2\nTTTT\n"             # plain header
    )
    seqs = bpc.parse_fasta_single(str(p))
    assert seqs["4471"] == "ACGTACGT"   # keyed by the locus, not the running index "1"
    assert seqs["locus2"] == "TTTT"


# ── parse_phase_output: empty-file contract (round-2 skip) ──────────────────
def test_parse_phase_output_empty_file(tmp_path):
    empty = tmp_path / "x.phase.out"
    empty.write_text("")
    n_chunks, _, _, _, _, _ = bpc.parse_phase_output(str(empty))
    assert n_chunks == 0


def test_parse_phase_output_missing_file(tmp_path):
    n_chunks, _, _, _, _, _ = bpc.parse_phase_output(str(tmp_path / "nope.phase.out"))
    assert n_chunks == 0


def test_parse_phase_output_one_block(tmp_path):
    p = tmp_path / "b.phase.out"
    p.write_text("#BLOCK: offset 1\n1\t0\t1\n2\t1\t0\n********\n")
    n_chunks, _, max_chunk, snp_matrix, chunks, ploidy = bpc.parse_phase_output(str(p))
    assert n_chunks == 1
    assert ploidy == 2
    assert chunks[1] == [1, 2]
    assert snp_matrix[1] == {1: "0", 2: "1"}
