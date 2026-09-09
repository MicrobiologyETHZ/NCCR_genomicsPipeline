"""
Unit tests for the PGAP samplesheet generator.

PGAP needs a per-genome taxon, which can't be inferred from a filename, so
this script writes an editable CSV rather than driving PGAP straight off a
directory scan. Pure Python — no conda environments, tools or databases are
needed.
"""
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow"))

from scripts.pgap_samplesheet import pgap_dir_to_samplesheet  # noqa: E402


def write_fasta(path, text=">contig_1\nACGT\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def read_csv_rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


@pytest.mark.unit
def test_writes_header_and_rows(tmp_path):
    d = tmp_path / "genomes"
    write_fasta(d / "strainA.fasta")
    out = tmp_path / "samples.csv"

    pgap_dir_to_samplesheet(d, out)

    rows = read_csv_rows(out)
    assert rows == [{"name": "strainA",
                     "fasta": str(d / "strainA.fasta"), "taxon": ""}]


@pytest.mark.unit
def test_scans_subdirectories(tmp_path):
    """The whole point vs. the phage assembly_dir default: nested genomes count."""
    d = tmp_path / "genomes"
    write_fasta(d / "top.fasta")
    write_fasta(d / "pacbio" / "nested.fasta")
    out = tmp_path / "samples.csv"

    pgap_dir_to_samplesheet(d, out)

    rows = read_csv_rows(out)
    assert {r["name"] for r in rows} == {"top", "nested"}


@pytest.mark.unit
def test_default_taxon_applied_to_every_row(tmp_path):
    d = tmp_path / "genomes"
    write_fasta(d / "strainA.fasta")
    write_fasta(d / "strainB.fasta")
    out = tmp_path / "samples.csv"

    pgap_dir_to_samplesheet(d, out, taxon="Pseudomonas")

    rows = read_csv_rows(out)
    assert all(r["taxon"] == "Pseudomonas" for r in rows)


@pytest.mark.unit
def test_gzipped_input_raises_clear_error(tmp_path):
    d = tmp_path / "genomes"
    (d / "compressed.fna.gz").parent.mkdir(parents=True, exist_ok=True)
    (d / "compressed.fna.gz").write_bytes(b"\x1f\x8b\x08\x00placeholder")
    out = tmp_path / "samples.csv"

    with pytest.raises(ValueError, match="decompress"):
        pgap_dir_to_samplesheet(d, out)


@pytest.mark.unit
def test_creates_output_directory(tmp_path):
    d = tmp_path / "genomes"
    write_fasta(d / "strainA.fasta")
    out = tmp_path / "nested" / "dir" / "samples.csv"

    pgap_dir_to_samplesheet(d, out)

    assert out.exists()


@pytest.mark.unit
def test_colliding_stems_disambiguated_by_parent_dir(tmp_path):
    """Reuses phage_inputs' collision handling rather than reinventing it."""
    d = tmp_path / "genomes"
    write_fasta(d / "runA" / "scaffolds.fasta")
    write_fasta(d / "runB" / "scaffolds.fasta")
    out = tmp_path / "samples.csv"

    pgap_dir_to_samplesheet(d, out)

    rows = read_csv_rows(out)
    assert {r["name"] for r in rows} == {"runA_scaffolds", "runB_scaffolds"}
