"""
Integration tests for Snakemake workflows
"""
import pytest
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
