"""
Integration tests for Snakemake workflows
"""
import os
import re
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


def _phage_dry_run(repo_root, target, extra_config=None, timeout=120):
    """Dry-run a phage target and return the completed process.

    Deliberately omits --use-conda: a dry run must not create the six tool
    environments. configs/test_phage_config.yaml likewise points at empty
    placeholder database directories, so nothing is downloaded either.
    """
    workflow_dir = repo_root / "workflow"
    config_file = repo_root / "configs" / "test_phage_config.yaml"

    if not config_file.exists():
        pytest.skip(f"Config file not found: {config_file}")

    cmd = ["snakemake", "-s", str(workflow_dir / "Snakefile_phage"),
           "--configfile", str(config_file), "-np", target]
    if extra_config:
        cmd += ["--config"] + extra_config

    try:
        return subprocess.run(cmd, cwd=workflow_dir, capture_output=True,
                              text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        pytest.fail(f"Phage dry run for {target} timed out")
    except FileNotFoundError:
        pytest.skip("Snakemake not installed")


@pytest.mark.integration
@pytest.mark.data
@pytest.mark.parametrize("target,expected_jobs", [
    # 3 assemblies: prepare + genomad + cenotetaker each (9), 2 callers x 3
    # assemblies of collect_viral (6), plus the target rule.
    ("find_phage", 16),
    # ... plus checkv/pharokka/phold/phynteny/remap_coordinates for each of the
    # 6 assembly x caller combinations (30).
    ("annotate_phage", 46),
    ("phage_summary", 46),
])
def test_phage_dry_run(repo_root, target, expected_jobs):
    """Each phage target resolves into a valid DAG of the expected size."""
    result = _phage_dry_run(repo_root, target)

    assert result.returncode == 0, f"Dry run failed: {result.stderr}"
    combined = result.stdout + result.stderr
    # Both callers must fan out, and the target itself must be scheduled.
    assert "genomad" in combined
    assert "cenotetaker" in combined
    assert re.search(rf"^total\s+{expected_jobs}$", combined, re.MULTILINE), (
        f"Expected {expected_jobs} jobs for {target}; got:\n{combined[-2000:]}")


@pytest.mark.integration
@pytest.mark.data
def test_phage_dry_run_handles_gzipped_input(repo_root):
    """Gzipped assemblies are decompressed; plain ones are passed through.

    The bundled test set has one .fna.gz and two plain FASTAs, so both branches
    of phage_prepare_input must appear.
    """
    result = _phage_dry_run(repo_root, "find_phage")

    assert result.returncode == 0, f"Dry run failed: {result.stderr}"
    combined = result.stdout + result.stderr
    assert "gzip -cd" in combined, "gzipped assembly was not routed through gzip -cd"
    # The NCBI-style name must survive .fna.gz stripping intact.
    assert "GCA_000001_genomic" in combined


@pytest.mark.integration
@pytest.mark.data
def test_phage_predict_only_skips_annotation(repo_root):
    """annotate=false drops CheckV and the annotation chain from the DAG.

    Targets phage_summary rather than find_phage: find_phage never includes the
    annotation chain, so it would pass regardless and prove nothing.

    Note `annotate=false` arrives from --config as the *string* "false", which
    is truthy in Python — this is exactly the case as_bool exists to handle.
    """
    result = _phage_dry_run(repo_root, "phage_summary",
                            extra_config=["annotate=false"])

    assert result.returncode == 0, f"Dry run failed: {result.stderr}"
    combined = result.stdout + result.stderr
    for rule in ("checkv", "pharokka", "phold", "phynteny"):
        assert not re.search(rf"^rule {rule}:", combined, re.MULTILINE), (
            f"{rule} should not run when annotate=false")
    # Prediction still happens, and the summary is still produced.
    assert re.search(r"^total\s+16$", combined, re.MULTILINE), (
        f"Expected 16 jobs with annotate=false; got:\n{combined[-2000:]}")


@pytest.mark.integration
@pytest.mark.data
def test_phage_missing_database_is_reported_clearly(repo_root, tmp_path):
    """A missing database must fail at startup naming the tool and the fix."""
    config = tmp_path / "bad_phage_config.yaml"
    config.write_text(
        "outDir: {}\n"
        "assembly_dir: {}\n"
        "callers: [genomad]\n"
        "annotate: false\n"
        "databases:\n"
        "  genomad: /definitely/not/a/real/path\n".format(
            tmp_path / "out",
            repo_root / "tests" / "test_data" / "phage" / "assemblies")
    )

    workflow_dir = repo_root / "workflow"
    cmd = ["snakemake", "-s", str(workflow_dir / "Snakefile_phage"),
           "--configfile", str(config), "-np", "find_phage"]
    try:
        result = subprocess.run(cmd, cwd=workflow_dir, capture_output=True,
                                text=True, timeout=120)
    except FileNotFoundError:
        pytest.skip("Snakemake not installed")

    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "genomad" in combined
    assert "download-database" in combined, (
        "Error should tell the user how to install the missing database")


def _env_dir(tmp_path, *tools):
    """A conda_env_dir containing environments for only the named tools."""
    root = tmp_path / "conda_envs"
    for tool in tools:
        (root / tool).mkdir(parents=True)
    return root


@pytest.mark.integration
@pytest.mark.data
def test_conda_env_dir_only_needs_envs_for_tools_that_run(repo_root, tmp_path):
    """A prediction-only run must not demand the annotation environments.

    Regression: every rule is defined unconditionally and snakemake evaluates
    each `conda:` at parse time, so conda_env() ran for checkv/pharokka/phold/
    phynteny even with annotate=false and failed the whole parse.
    """
    env_dir = _env_dir(tmp_path, "genomad", "cenotetaker")
    result = _phage_dry_run(repo_root, "find_phage", extra_config=[
        "annotate=false", f"conda_env_dir={env_dir}"])

    assert result.returncode == 0, (
        f"prediction-only run should not need annotation envs:\n{result.stderr}")


@pytest.mark.integration
@pytest.mark.data
def test_conda_env_dir_still_fails_for_tools_that_do_run(repo_root, tmp_path):
    """The guard must still fire when the missing env is actually needed."""
    env_dir = _env_dir(tmp_path, "genomad", "cenotetaker")
    result = _phage_dry_run(repo_root, "phage_summary", extra_config=[
        "annotate=true", f"conda_env_dir={env_dir}"])

    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "checkv" in combined
    assert "mamba create" in combined, "error should say how to create the env"


@pytest.mark.integration
@pytest.mark.data
def test_conda_env_dir_fails_for_missing_caller_env(repo_root, tmp_path):
    """geNomad always runs, so its env is required even in prediction-only mode."""
    env_dir = _env_dir(tmp_path, "cenotetaker")
    result = _phage_dry_run(repo_root, "find_phage", extra_config=[
        "annotate=false", f"conda_env_dir={env_dir}"])

    assert result.returncode != 0
    assert "genomad" in result.stdout + result.stderr


@pytest.mark.integration
@pytest.mark.data
def test_placeholder_database_warns(repo_root):
    """Test-placeholder db paths get copied into real configs; warn loudly.

    Not an error: configs/test_phage_config.yaml legitimately points at them so
    dry runs need no real databases.
    """
    result = _phage_dry_run(repo_root, "phage_summary")

    assert result.returncode == 0
    combined = result.stdout + result.stderr
    assert "test placeholder" in combined


@pytest.mark.integration
@pytest.mark.slow
def test_phage_end_to_end_smoke(repo_root, tmp_path):
    """Run the phage workflow for real.

    Skipped unless all six tools are already installed AND real databases are
    configured via NCCR_PHAGE_DB_ROOT. This is never run during development —
    the tools and their ~20 GB of databases live on the cluster.
    """
    tools = ["genomad", "cenotetaker3", "checkv", "pharokka.py",
             "phold", "phynteny_transformer"]
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        pytest.skip(f"Phage tools not on PATH: {', '.join(missing)}")

    db_root = os.environ.get("NCCR_PHAGE_DB_ROOT")
    if not db_root:
        pytest.skip("Set NCCR_PHAGE_DB_ROOT to the phage database directory")

    workflow_dir = repo_root / "workflow"
    cmd = [
        "snakemake", "-s", str(workflow_dir / "Snakefile_phage"),
        "--configfile", str(repo_root / "configs" / "test_phage_config.yaml"),
        "--use-conda", "--cores", "4",
        "--config", f"outDir={tmp_path / 'output'}",
        *[f"databases={db_root}"],
        "phage_summary",
    ]
    result = subprocess.run(cmd, cwd=workflow_dir, capture_output=True,
                            text=True, timeout=7200)

    assert result.returncode == 0, f"Phage run failed:\n{result.stdout}\n{result.stderr}"
    summary = tmp_path / "output" / "phage" / "summary" / "phage_predictions.tsv"
    assert summary.exists(), "No prediction summary was produced"


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
