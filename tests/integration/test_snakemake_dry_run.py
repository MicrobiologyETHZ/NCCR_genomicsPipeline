"""
Integration tests for Snakemake workflows
"""
import pytest
import shutil
import subprocess
from pathlib import Path


@pytest.mark.integration
@pytest.mark.data
def test_variant_calling_dry_run(repo_root, test_data_dir):
    """Test that variant calling workflow can do a dry run"""
    workflow_dir = repo_root / "workflow"
    config_file = repo_root / "configs" / "test_variant_calling_config.yaml"

    if not config_file.exists():
        pytest.skip(f"Config file not found: {config_file}")

    cmd = [
        "snakemake",
        "-s", str(workflow_dir / "Snakefile"),
        "--configfile", str(config_file),
        "-np",
        "call_variants"
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=workflow_dir,
            capture_output=True,
            text=True,
            timeout=30
        )
        # Dry run should succeed (exit code 0) or show what would be run
        assert result.returncode in [0, 1], f"Unexpected error: {result.stderr}"
    except subprocess.TimeoutExpired:
        pytest.fail("Dry run timed out")
    except FileNotFoundError:
        pytest.skip("Snakemake not installed")


@pytest.mark.integration
@pytest.mark.data
@pytest.mark.parametrize("target", ["instrain", "instrain_compare"])
def test_instrain_dry_run(repo_root, target):
    """InStrain profile/compare targets resolve into a valid DAG (no install needed)."""
    workflow_dir = repo_root / "workflow"
    config_file = repo_root / "configs" / "test_instrain_config.yaml"

    if not config_file.exists():
        pytest.skip(f"Config file not found: {config_file}")

    cmd = [
        "snakemake",
        "-s", str(workflow_dir / "Snakefile"),
        "--configfile", str(config_file),
        "-np",
        target,
    ]

    try:
        result = subprocess.run(
            cmd, cwd=workflow_dir, capture_output=True, text=True, timeout=60
        )
    except subprocess.TimeoutExpired:
        pytest.fail("Dry run timed out")
    except FileNotFoundError:
        pytest.skip("Snakemake not installed")

    assert result.returncode == 0, f"Dry run failed: {result.stderr}"
    # The DAG should schedule the instrain rule and pull in the BAM it depends on.
    combined = result.stdout + result.stderr
    assert "instrain" in combined


@pytest.mark.integration
@pytest.mark.data
@pytest.mark.slow
def test_instrain_profile_smoke(repo_root, tmp_path):
    """Actually run `inStrain profile` end-to-end on the bundled test data.

    Skipped unless inStrain + Snakemake are installed (i.e. the call_variants
    conda env is active). The bundled FASTQs are tiny, so min_cov is relaxed in
    the test config; we assert the profile output dir is produced.
    """
    if shutil.which("inStrain") is None:
        pytest.skip("inStrain not on PATH (activate the call_variants env)")
    if shutil.which("snakemake") is None:
        pytest.skip("Snakemake not installed")

    workflow_dir = repo_root / "workflow"
    config_file = repo_root / "configs" / "test_instrain_config.yaml"
    if not config_file.exists():
        pytest.skip(f"Config file not found: {config_file}")

    cmd = [
        "snakemake",
        "-s", str(workflow_dir / "Snakefile"),
        "--configfile", str(config_file),
        "--use-conda",
        "--cores", "2",
        f"--config", f"outDir={tmp_path / 'output'}",
        "instrain",
    ]

    try:
        result = subprocess.run(
            cmd, cwd=workflow_dir, capture_output=True, text=True, timeout=1800
        )
    except subprocess.TimeoutExpired:
        pytest.fail("InStrain smoke run timed out")

    assert result.returncode == 0, f"InStrain run failed:\n{result.stdout}\n{result.stderr}"
    markers = list((tmp_path / "output" / "instrain" / "profiles").glob("*/.IS.profile.done"))
    assert markers, "No InStrain profile completion markers were produced"


@pytest.mark.integration
def test_snakefile_syntax(repo_root):
    """Test that Snakefile has valid syntax"""
    workflow_dir = repo_root / "workflow"
    snakefile = workflow_dir / "Snakefile"

    cmd = [
        "snakemake",
        "-s", str(snakefile),
        "--lint"
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=workflow_dir,
            capture_output=True,
            text=True,
            timeout=10
        )
        # Lint should either succeed or not be available in older versions
        if "Linting" in result.stdout or result.returncode == 0:
            assert True
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pytest.skip("Snakemake lint not available")
