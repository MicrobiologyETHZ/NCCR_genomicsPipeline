"""
Unit tests for utility functions
"""
import pytest
from pathlib import Path


@pytest.mark.unit
def test_repo_structure(repo_root):
    """Test that expected directories exist"""
    assert (repo_root / "workflow").exists()
    assert (repo_root / "workflow" / "Snakefile").exists()
    assert (repo_root / "workflow" / "main.py").exists()
    assert (repo_root / "workflow" / "rules").exists()
    assert (repo_root / "workflow" / "envs").exists()


@pytest.mark.unit
def test_test_data_exists(test_data_dir):
    """Test that test data is present"""
    assert test_data_dir.exists()
    assert (test_data_dir / "test_samples.txt").exists()
    assert (test_data_dir / "LL6_1.fasta.gz").exists()


@pytest.mark.unit
def test_sample_file_readable(sample_list):
    """Test that sample file can be read"""
    assert len(sample_list) > 0
    assert "LL23" in sample_list
