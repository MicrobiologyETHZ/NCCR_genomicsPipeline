"""
Unit tests for mapping phage annotations back to genome coordinates.

Getting this wrong is a silent failure: an off-by-one or a reverse-complemented
extract yields coordinates that look entirely plausible but point at the wrong
place. The fixtures below use real values from a Leaf257 run, where the
extraction offsets were confirmed independently with `samtools faidx`.

Pure data transformation — no tools, environments or databases needed.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow"))

from scripts.phage_collect import (  # noqa: E402
    cenotetaker_coords,
    genomad_coords,
)
from scripts.phage_remap_coords import (  # noqa: E402
    parse_attributes,
    parse_gff,
    remap,
    verify_offsets,
)


# Real rows from the Leaf257 run.
GENOMAD_SUMMARY = (
    "seq_name\tlength\ttopology\tcoordinates\tn_genes\n"
    "contig_153|provirus_907014_947422\t40409\tProvirus\t907014-947422\t43\n"
    "contig_153|provirus_1_16548\t16548\tProvirus\t1-16548\t20\n"
    "contig_99\t45000\tNo terminal repeats\tNA\t50\n"
)
CT_VIRUS_SUMMARY = (
    "contig\tinput_name\tvirus_seq_length\n"
    "Leaf257_5@C38\tcontig_153\t39618\n"
)
CT_PRUNE_SUMMARY = (
    "contig\tcontig_length\tchunk_length\tchunk_name\tchunk_start\tchunk_stop\n"
    "Leaf257_5\t3964747\t39618\tC38\t710296\t749914\n"
)


def write(path, text):
    path.write_text(text)
    return str(path)


@pytest.fixture
def genomad(tmp_path):
    return genomad_coords(write(tmp_path / "gen.tsv", GENOMAD_SUMMARY))


@pytest.fixture
def cenotetaker(tmp_path):
    return cenotetaker_coords(
        write(tmp_path / "virus.tsv", CT_VIRUS_SUMMARY),
        write(tmp_path / "prune.tsv", CT_PRUNE_SUMMARY),
    )


# --------------------------------------------------------------- caller parsing

@pytest.mark.unit
def test_genomad_coordinates_are_1indexed_inclusive(genomad):
    """Verified by geNomad's own length column: end - start + 1 == length."""
    row = genomad["contig_153|provirus_907014_947422"]

    assert row["source_contig"] == "contig_153"
    assert row["source_start"] == 907014
    assert row["source_end"] == 947422
    assert row["offset"] == 907013          # 1-indexed start -> 0-based offset
    assert row["source_end"] - row["source_start"] + 1 == 40409


@pytest.mark.unit
def test_genomad_provirus_starting_at_base_one_has_zero_offset(genomad):
    assert genomad["contig_153|provirus_1_16548"]["offset"] == 0


@pytest.mark.unit
def test_genomad_whole_contig_virus_has_no_offset(genomad):
    """`coordinates` is NA when the virus is not integrated."""
    row = genomad["contig_99"]

    assert row["offset"] == 0
    assert row["is_provirus"] is False
    assert row["source_contig"] == "contig_99"


@pytest.mark.unit
def test_cenotetaker_chunk_start_is_0indexed_half_open(cenotetaker):
    """Different convention from geNomad; confirmed with samtools faidx.

    chunk_stop - chunk_start == chunk_length (not +1), and the extract's first
    base sits at 1-based 710297, i.e. chunk_start + 1.
    """
    row = cenotetaker["Leaf257_5@C38"]

    assert row["source_contig"] == "contig_153"   # via input_name, not the ID
    assert row["offset"] == 710296
    assert row["offset"] + 1 == 710297
    assert 749914 - 710296 == 39618


@pytest.mark.unit
def test_both_callers_agree_on_the_same_prophage_start(genomad, tmp_path):
    """geNomad and Cenote-Taker placed the same prophage at the same base.

    geNomad: provirus_710297_749031 (1-based). CT3: chunk_start 710296 (0-based).
    Both must normalise to the same genomic position for the first base.
    """
    gen = genomad_coords(write(tmp_path / "g.tsv",
        "seq_name\tlength\ttopology\tcoordinates\tn_genes\n"
        "contig_153|provirus_710297_749031\t38735\tProvirus\t710297-749031\t40\n"))
    ct = cenotetaker_coords(
        write(tmp_path / "v.tsv", CT_VIRUS_SUMMARY),
        write(tmp_path / "p.tsv", CT_PRUNE_SUMMARY))

    first_base_genomad = gen["contig_153|provirus_710297_749031"]["offset"] + 1
    first_base_ct = ct["Leaf257_5@C38"]["offset"] + 1

    assert first_base_genomad == first_base_ct == 710297


# ------------------------------------------------------------------- transform

def _coords_row(**overrides):
    row = {
        "viral_id": "v1", "assembly": "Leaf257", "caller": "genomad",
        "source_contig": "contig_153", "offset": "907013",
        "source_start": "907014", "source_end": "947422",
        "extracted_length": "40409", "is_provirus": "True", "orientation": "+",
    }
    row.update(overrides)
    return {row["viral_id"]: row}


def _gff(tmp_path, *feature_lines):
    path = tmp_path / "pharokka.gff"
    path.write_text("##gff-version 3\n" + "".join(feature_lines))
    return str(path)


@pytest.mark.unit
def test_worked_example(tmp_path):
    """A CDS at 1200-1800 inside a provirus at 907014 lands at 908213."""
    gff = _gff(tmp_path,
               "v1\tPyrodigal\tCDS\t1200\t1800\t0\t+\t0\tID=v1_CDS_0001;product=tail\n")

    rows, unmapped = remap(gff, _coords_row(), "Leaf257", "genomad")

    assert not unmapped
    assert rows[0]["genome_start"] == 908213
    assert rows[0]["genome_end"] == 908813
    assert rows[0]["source_contig"] == "contig_153"


@pytest.mark.unit
def test_first_base_of_extract_maps_to_source_start(tmp_path):
    gff = _gff(tmp_path, "v1\tPyrodigal\tCDS\t1\t100\t0\t+\t0\tID=a\n")

    rows, _ = remap(gff, _coords_row(), "Leaf257", "genomad")

    assert rows[0]["genome_start"] == 907014


@pytest.mark.unit
def test_last_base_of_extract_maps_to_source_end(tmp_path):
    """40409 is the full extract length, so its last base is source_end."""
    gff = _gff(tmp_path, "v1\tPyrodigal\tCDS\t40409\t40409\t0\t+\t0\tID=a\n")

    rows, _ = remap(gff, _coords_row(), "Leaf257", "genomad")

    assert rows[0]["genome_end"] == 947422


@pytest.mark.unit
def test_zero_offset_leaves_coordinates_unchanged(tmp_path):
    gff = _gff(tmp_path, "v1\tPyrodigal\tCDS\t500\t900\t0\t+\t0\tID=a\n")

    rows, _ = remap(gff, _coords_row(offset="0", source_start="1"),
                    "Leaf257", "genomad")

    assert (rows[0]["genome_start"], rows[0]["genome_end"]) == (500, 900)


@pytest.mark.unit
def test_unmapped_sequence_is_reported_not_silently_dropped(tmp_path):
    gff = _gff(tmp_path, "ghost\tPyrodigal\tCDS\t1\t100\t0\t+\t0\tID=a\n")

    rows, unmapped = remap(gff, _coords_row(), "Leaf257", "genomad")

    assert rows == []
    assert unmapped == ["ghost"]


@pytest.mark.unit
def test_legacy_assembly_prefix_on_seqid_still_resolves(tmp_path):
    """Older runs prefixed FASTA headers with '<assembly>|'."""
    gff = _gff(tmp_path, "Leaf257|v1\tPyrodigal\tCDS\t1\t100\t0\t+\t0\tID=a\n")

    rows, unmapped = remap(gff, _coords_row(), "Leaf257", "genomad")

    assert not unmapped
    assert rows[0]["genome_start"] == 907014


@pytest.mark.unit
def test_non_cds_features_are_carried_through(tmp_path):
    """pharokka emits tRNA and CRISPR features too; all need remapping."""
    gff = _gff(tmp_path,
               "v1\tPyrodigal\tCDS\t1\t100\t0\t+\t0\tID=a\n"
               "v1\ttRNAscan-SE\ttRNA\t200\t280\t0\t-\t.\tID=b\n")

    rows, _ = remap(gff, _coords_row(), "Leaf257", "genomad")

    assert [r["feature_type"] for r in rows] == ["CDS", "tRNA"]


@pytest.mark.unit
def test_tab_inside_attributes_does_not_break_parsing(tmp_path):
    """Observed in pharokka tRNA rows: a raw tab inside the attributes column.

    A naive split produces ten fields and mangles the record.
    """
    gff = _gff(tmp_path,
               "v1\ttRNAscan-SE\ttRNA\t1965\t2046\t69.1\t-\t.\t"
               "ID=v1_tRNA_0001;product=tRNA-Leu(UUG);"
               "anticodon=(pos:2012..2010,aa:Leu,seq:CAA)\tUUG\n")

    features = list(parse_gff(gff))

    assert len(features) == 1
    assert features[0]["start"] == 1965
    assert parse_attributes(features[0]["attributes"])["product"] == "tRNA-Leu(UUG)"


# ---------------------------------------------------------------- verification

@pytest.mark.unit
def test_verify_offsets_accepts_a_correct_offset(tmp_path):
    contig = "AAAACCCCGGGGTTTT" * 10
    (tmp_path / "asm.fna").write_text(f">contig_153\n{contig}\n")
    (tmp_path / "viral.fna").write_text(f">v1\n{contig[100:150]}\n")

    problems = verify_offsets(_coords_row(offset="100", source_contig="contig_153"),
                              str(tmp_path / "viral.fna"),
                              str(tmp_path / "asm.fna"))

    assert problems == []


@pytest.mark.unit
def test_verify_offsets_catches_off_by_one(tmp_path):
    """The exact failure this guard exists for: coordinates that look fine."""
    contig = "".join("ACGT"[i % 4] for i in range(400))
    contig = contig[:100] + "TTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAA" + contig[140:]
    (tmp_path / "asm.fna").write_text(f">contig_153\n{contig}\n")
    (tmp_path / "viral.fna").write_text(f">v1\n{contig[100:150]}\n")

    problems = verify_offsets(_coords_row(offset="101", source_contig="contig_153"),
                              str(tmp_path / "viral.fna"),
                              str(tmp_path / "asm.fna"))

    assert problems and "does not match" in problems[0]


@pytest.mark.unit
def test_verify_offsets_catches_reverse_complement(tmp_path):
    """Cenote-Taker's pipeline reorients contigs; a flipped extract must fail."""
    contig = "".join("ACGT"[i % 4] for i in range(200))
    region = contig[50:100]
    revcomp = region[::-1].translate(str.maketrans("ACGT", "TGCA"))
    (tmp_path / "asm.fna").write_text(f">contig_153\n{contig}\n")
    (tmp_path / "viral.fna").write_text(f">v1\n{revcomp}\n")

    problems = verify_offsets(_coords_row(offset="50", source_contig="contig_153"),
                              str(tmp_path / "viral.fna"),
                              str(tmp_path / "asm.fna"))

    assert problems, "a reverse-complemented extract must not pass verification"


@pytest.mark.unit
def test_verify_offsets_is_a_noop_without_fastas(tmp_path):
    assert verify_offsets(_coords_row(), str(tmp_path / "absent.fna"),
                          str(tmp_path / "absent2.fna")) == []
