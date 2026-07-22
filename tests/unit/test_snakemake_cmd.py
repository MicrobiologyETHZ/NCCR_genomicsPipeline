"""
Unit tests for the snakemake command the CLI builds.

Snakemake 8 deprecated and 9 removed `--cluster`; submission now goes through an
executor plugin. These tests pin the resulting invocation so the regression that
broke every cluster command ("unrecognized arguments: --cluster") cannot recur.

Nothing is executed: subprocess.check_call is intercepted, so these run without
the executor plugin installed, without SLURM, and without any tool environment.
"""
import sys
from pathlib import Path
from unittest import mock

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from workflow import main as cli  # noqa: E402


CONFIG = "configs/test_phage_config.yaml"


def build(**kwargs):
    """Return the argv list snakemake_cmd would run, without running it."""
    kwargs.setdefault("config", CONFIG)
    kwargs.setdefault("analysis", "find_phage")
    kwargs.setdefault("smk_file", "Snakefile_phage")
    kwargs.setdefault("dry", False)
    kwargs.setdefault("local", False)
    with mock.patch.object(cli.subprocess, "check_call") as call:
        # Pretend the executor plugin is installed; its absence is tested below.
        with mock.patch.object(cli, "check_cluster_executor"):
            cli.snakemake_cmd(**kwargs)
    return [str(a) for a in call.call_args[0][0]]


@pytest.mark.unit
def test_cluster_run_uses_executor_plugin_not_removed_cluster_flag():
    argv = build()

    assert "--executor" in argv
    assert argv[argv.index("--executor") + 1] == "cluster-generic"
    assert "--cluster-generic-submit-cmd" in argv
    # The flag removed in Snakemake 9 must not reappear.
    assert "--cluster" not in argv


@pytest.mark.unit
def test_submit_command_is_a_single_argument():
    """The sbatch string must reach snakemake as one argv element.

    Splitting it would make snakemake read `sbatch`, `-t`, ... as separate
    options and fail.
    """
    argv = build()
    submit = argv[argv.index("--cluster-generic-submit-cmd") + 1]

    assert submit.startswith("DIR=$(dirname")
    assert "sbatch" in submit


@pytest.mark.unit
def test_submit_command_keeps_snakemake_placeholders_unexpanded():
    """Resources come from each rule's params; snakemake expands them per job."""
    argv = build(partition="institute")
    submit = argv[argv.index("--cluster-generic-submit-cmd") + 1]

    for placeholder in ("{params.time}", "{params.mem}", "{threads}",
                        "{params.qoutfile}", "{params.qerrfile}"):
        assert placeholder in submit, f"{placeholder} was lost or expanded"
    assert "--partition institute" in submit


@pytest.mark.unit
def test_submit_command_has_no_stray_braces():
    """Any brace that is not a placeholder would break snakemake's formatting."""
    import re

    submit = build()[0:]
    submit = submit[submit.index("--cluster-generic-submit-cmd") + 1]
    # Strip the known placeholders, then assert nothing brace-like remains.
    stripped = re.sub(r"\{(params\.\w+|threads)\}", "", submit)
    assert "{" not in stripped and "}" not in stripped


@pytest.mark.unit
def test_partition_is_configurable():
    argv = build(partition="gpu")
    submit = argv[argv.index("--cluster-generic-submit-cmd") + 1]
    assert "--partition gpu" in submit


@pytest.mark.unit
def test_latency_wait_passed_for_shared_filesystem():
    """Output lands on NFS, where a job can finish before its files appear."""
    argv = build(latency_wait=120)
    assert argv[argv.index("--latency-wait") + 1] == "120"


@pytest.mark.unit
@pytest.mark.parametrize("no_conda,expected", [(False, True), (True, False)])
def test_use_conda_honours_no_conda_on_cluster(no_conda, expected):
    argv = build(no_conda=no_conda)
    assert ("--use-conda" in argv) is expected


@pytest.mark.unit
@pytest.mark.parametrize("no_conda,expected", [(False, True), (True, False)])
def test_use_conda_honours_no_conda_locally(no_conda, expected):
    """--local previously ignored conda entirely; both paths must agree."""
    argv = build(local=True, no_conda=no_conda)
    assert ("--use-conda" in argv) is expected


@pytest.mark.unit
def test_cores_default_and_override():
    assert build()[build().index("-j") + 1] == "6"          # cluster default
    assert build(cores=20)[build(cores=20).index("-j") + 1] == "20"
    local = build(local=True)
    assert local[local.index("-j") + 1] == "1"              # local default


@pytest.mark.unit
def test_dry_run_is_unchanged_and_has_no_cluster_flags():
    argv = build(dry=True)

    assert "-np" in argv
    assert "--executor" not in argv
    assert "--cluster-generic-submit-cmd" not in argv


@pytest.mark.unit
def test_config_overrides_are_appended_last():
    argv = build(set_config=("annotate=false", "callers=genomad"))

    assert argv[-3:] == ["--config", "annotate=false", "callers=genomad"]


@pytest.mark.unit
@pytest.mark.parametrize("analysis,smk", [
    ("assemble", "Snakefile"),
    ("breseq", "Snakefile"),
    ("instrain", "Snakefile"),
    ("ismap", "Snakefile"),
    ("phage_summary", "Snakefile_phage"),
])
def test_all_commands_share_the_fixed_cluster_path(analysis, smk):
    """The break was in shared code, so verify beyond the command that hit it."""
    argv = build(analysis=analysis, smk_file=smk)

    assert "--cluster" not in argv
    assert "--cluster-generic-submit-cmd" in argv
    assert argv[-1] == analysis


@pytest.mark.unit
def test_missing_executor_plugin_gives_actionable_error():
    """The real failure mode: snakemake's own message names no remedy."""
    with mock.patch.object(cli.importlib.util, "find_spec", return_value=None):
        with pytest.raises(Exception) as excinfo:
            cli.check_cluster_executor()

    message = str(excinfo.value)
    assert "snakemake-executor-plugin-cluster-generic" in message
    assert "conda install" in message


@pytest.mark.unit
def test_executor_check_passes_when_plugin_present():
    with mock.patch.object(cli.importlib.util, "find_spec", return_value=object()):
        cli.check_cluster_executor()


@pytest.mark.unit
@pytest.mark.parametrize("kwargs", [{"dry": True}, {"local": True}])
def test_executor_check_not_triggered_for_dry_or_local(kwargs):
    """Neither path submits jobs, so neither needs the plugin installed."""
    with mock.patch.object(cli.subprocess, "check_call"):
        with mock.patch.object(cli, "check_cluster_executor") as check:
            cli.snakemake_cmd(CONFIG, "find_phage", "Snakefile_phage",
                              dry=kwargs.get("dry", False),
                              local=kwargs.get("local", False))
    check.assert_not_called()


@pytest.mark.unit
def test_missing_config_file_raises_before_anything_runs():
    with mock.patch.object(cli.subprocess, "check_call") as call:
        with pytest.raises(FileNotFoundError):
            cli.snakemake_cmd("configs/does_not_exist.yaml", "find_phage",
                              "Snakefile_phage", dry=True, local=False)
    call.assert_not_called()
