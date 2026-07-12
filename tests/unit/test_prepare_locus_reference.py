"""Unit tests for bin/prepare_locus_reference.py."""
import prepare_locus_reference as plr


def test_parse_fasta_keys_by_header(tmp_path):
    p = tmp_path / "locus1.fasta"
    p.write_text(">sample1\nACGT\n>sample2\nTTTT\n")
    seqs = plr.parse_fasta(str(p))
    assert seqs == {"sample1": "ACGT", "sample2": "TTTT"}


def test_parse_fasta_strips_hybpiper_suffix(tmp_path):
    p = tmp_path / "locus1.fasta"
    p.write_text(">sample1-L001\nACGTACGT\n")
    seqs = plr.parse_fasta(str(p))
    assert "sample1" in seqs               # -L001 suffix stripped
    assert seqs["sample1"] == "ACGTACGT"


def test_parse_fasta_multiline_sequence(tmp_path):
    p = tmp_path / "locus1.fasta"
    p.write_text(">sample1\nACGT\nACGT\nAC\n")
    seqs = plr.parse_fasta(str(p))
    assert seqs["sample1"] == "ACGTACGTAC"
