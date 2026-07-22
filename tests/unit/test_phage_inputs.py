"""
Unit tests for phage workflow input resolution.

These cover the logic that maps arbitrary assembly files to unique short names,
which is what lets the phage workflow run on public/collaborator data rather
than only on assemblies produced by this pipeline. Pure Python — no conda
environments, tools or databases are needed.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow"))

from scripts.phage_inputs import (  # noqa: E402
    as_bool,
    assembly_name,
    is_compressed,
    resolve_assemblies,
    validate_databases,
)


@pytest.fixture
def basedir(tmp_path):
    """Stand-in for the workflow/ directory that relative paths resolve against."""
    return tmp_path


def write_fasta(path, text=">contig_1\nACGT\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.mark.unit
@pytest.mark.parametrize("filename,expected", [
    ("scaffolds.fasta", "scaffolds"),
    ("strainA.fa", "strainA"),
    ("GCA_000001_genomic.fna", "GCA_000001_genomic"),
    ("GCA_000001_genomic.fna.gz", "GCA_000001_genomic"),
    ("assembly.fasta.gz", "assembly"),
    # A dotted name that isn't a known suffix must survive intact.
    ("Ecoli.K12.fna", "Ecoli.K12"),
])
def test_assembly_name_strips_suffixes(filename, expected):
    assert assembly_name(filename) == expected


@pytest.mark.unit
@pytest.mark.parametrize("value,expected", [
    (True, True),
    (False, False),
    # Snakemake's `--config annotate=false` yields the STRING "false", which is
    # truthy in Python. Getting this wrong silently runs the whole annotation
    # chain when the user asked for prediction only.
    ("false", False),
    ("False", False),
    ("no", False),
    ("off", False),
    ("0", False),
    ("", False),
    ("true", True),
    ("yes", True),
])
def test_as_bool_handles_config_strings(value, expected):
    assert as_bool(value) is expected


@pytest.mark.unit
def test_as_bool_default_used_only_for_none():
    assert as_bool(None, default=True) is True
    assert as_bool(None, default=False) is False
    assert as_bool("false", default=True) is False


@pytest.mark.unit
def test_is_compressed():
    assert is_compressed("a.fna.gz")
    assert not is_compressed("a.fna")


@pytest.mark.unit
def test_resolve_from_assembly_dir(tmp_path, basedir):
    d = tmp_path / "assemblies"
    write_fasta(d / "strainA.fasta")
    write_fasta(d / "strainB.fna")
    (d / "notes.txt").write_text("ignore me")

    got = resolve_assemblies({"assembly_dir": str(d)}, basedir)

    assert set(got) == {"strainA", "strainB"}
    assert got["strainA"] == d / "strainA.fasta"


@pytest.mark.unit
def test_assembly_dir_matches_gzipped_without_explicit_pattern(tmp_path, basedir):
    """Public assemblies usually arrive gzipped; the default patterns must catch them."""
    d = tmp_path / "assemblies"
    write_fasta(d / "plain.fna")
    (d / "compressed.fna.gz").write_bytes(b"\x1f\x8b\x08\x00placeholder")

    got = resolve_assemblies({"assembly_dir": str(d)}, basedir)

    assert set(got) == {"plain", "compressed"}


@pytest.mark.unit
def test_explicit_pattern_narrows_selection(tmp_path, basedir):
    d = tmp_path / "assemblies"
    write_fasta(d / "keep.fna")
    write_fasta(d / "skip.fasta")

    got = resolve_assemblies({"assembly_dir": str(d), "pattern": "*.fna"}, basedir)

    assert set(got) == {"keep"}


@pytest.mark.unit
def test_resolve_from_assemblies_list(tmp_path, basedir):
    a = write_fasta(tmp_path / "one" / "asmA.fna")
    b = write_fasta(tmp_path / "two" / "asmB.fna")

    got = resolve_assemblies({"assemblies": [str(a), str(b)]}, basedir)

    assert got == {"asmA": a, "asmB": b}


@pytest.mark.unit
def test_resolve_from_assemblies_mapping_controls_names(tmp_path, basedir):
    a = write_fasta(tmp_path / "runA" / "scaffolds.fasta")
    b = write_fasta(tmp_path / "runB" / "scaffolds.fasta")

    got = resolve_assemblies(
        {"assemblies": {"strainA": str(a), "strainB": str(b)}}, basedir)

    assert got == {"strainA": a, "strainB": b}


@pytest.mark.unit
def test_colliding_stems_fall_back_to_parent_dir(tmp_path, basedir):
    """Several pipeline runs all produce `scaffolds.fasta` — names must not clash."""
    a = write_fasta(tmp_path / "runA" / "scaffolds.fasta")
    b = write_fasta(tmp_path / "runB" / "scaffolds.fasta")

    got = resolve_assemblies({"assemblies": [str(a), str(b)]}, basedir)

    assert set(got) == {"runA_scaffolds", "runB_scaffolds"}
    assert got["runA_scaffolds"] == a
    assert got["runB_scaffolds"] == b


@pytest.mark.unit
def test_unresolvable_collision_raises(tmp_path, basedir):
    """Same filename AND same parent dir name: must fail loudly, not overwrite."""
    a = write_fasta(tmp_path / "x" / "run" / "scaffolds.fasta")
    b = write_fasta(tmp_path / "y" / "run" / "scaffolds.fasta")

    with pytest.raises(ValueError, match="unique names"):
        resolve_assemblies({"assemblies": [str(a), str(b)]}, basedir)


@pytest.mark.unit
def test_relative_paths_resolve_against_basedir(tmp_path):
    write_fasta(tmp_path / "assemblies" / "rel.fna")

    got = resolve_assemblies({"assembly_dir": "assemblies"}, tmp_path)

    assert got["rel"] == tmp_path / "assemblies" / "rel.fna"


@pytest.mark.unit
def test_sources_combine(tmp_path, basedir):
    d = tmp_path / "assemblies"
    write_fasta(d / "fromdir.fna")
    extra = write_fasta(tmp_path / "extra" / "explicit.fna")

    got = resolve_assemblies(
        {"assembly_dir": str(d), "assemblies": [str(extra)]}, basedir)

    assert set(got) == {"fromdir", "explicit"}


@pytest.mark.unit
def test_no_sources_raises(basedir):
    with pytest.raises(ValueError, match="No assemblies to process"):
        resolve_assemblies({}, basedir)


@pytest.mark.unit
def test_missing_assembly_dir_raises(tmp_path, basedir):
    with pytest.raises(ValueError, match="does not exist"):
        resolve_assemblies({"assembly_dir": str(tmp_path / "nope")}, basedir)


@pytest.mark.unit
def test_empty_assembly_dir_raises(tmp_path, basedir):
    (tmp_path / "empty").mkdir()
    with pytest.raises(ValueError, match="No assemblies found"):
        resolve_assemblies({"assembly_dir": str(tmp_path / "empty")}, basedir)


@pytest.mark.unit
def test_validate_databases_accepts_existing_paths(tmp_path, basedir):
    db = tmp_path / "genomad_db"
    db.mkdir()

    got = validate_databases(
        {"databases": {"genomad": str(db)}}, {"genomad"}, basedir)

    assert got == {"genomad": str(db)}


@pytest.mark.unit
def test_validate_databases_reports_missing_path_and_install_command(basedir):
    with pytest.raises(ValueError) as excinfo:
        validate_databases({"databases": {}}, {"phold"}, basedir)

    message = str(excinfo.value)
    assert "phold" in message
    # The error must tell the user how to fix it, not just that it's broken.
    assert "phold install" in message


@pytest.mark.unit
def test_validate_databases_rejects_nonexistent_directory(tmp_path, basedir):
    with pytest.raises(ValueError, match="does not exist"):
        validate_databases(
            {"databases": {"checkv": str(tmp_path / "absent")}}, {"checkv"}, basedir)


@pytest.mark.unit
def test_validate_databases_warns_on_test_placeholder(tmp_path, basedir, capsys):
    """Placeholder paths copied out of the test config must be called out.

    A warning rather than an error: configs/test_phage_config.yaml points at
    these on purpose so dry runs need no real databases.
    """
    db = tmp_path / "genomad_db"
    db.mkdir()
    (db / ".placeholder").write_text("not a real database")

    got = validate_databases(
        {"databases": {"genomad": str(db)}}, {"genomad"}, basedir)

    assert got == {"genomad": str(db)}          # still usable, just flagged
    warning = capsys.readouterr().err
    assert "test placeholder" in warning
    assert "genomad download-database" in warning


@pytest.mark.unit
def test_validate_databases_silent_for_real_looking_database(tmp_path, basedir, capsys):
    db = tmp_path / "genomad_db"
    db.mkdir()
    (db / "genomad_db.dmnd").write_text("pretend index")

    validate_databases({"databases": {"genomad": str(db)}}, {"genomad"}, basedir)

    assert "placeholder" not in capsys.readouterr().err


@pytest.mark.unit
def test_validate_databases_only_checks_required_tools(tmp_path, basedir):
    """A prediction-only run must not demand the annotation databases."""
    db = tmp_path / "genomad_db"
    db.mkdir()

    got = validate_databases(
        {"databases": {"genomad": str(db)}}, {"genomad"}, basedir)

    assert set(got) == {"genomad"}
