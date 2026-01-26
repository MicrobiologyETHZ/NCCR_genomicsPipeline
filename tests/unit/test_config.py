"""
Unit tests for configuration handling
"""
import pytest
import yaml
from pathlib import Path


@pytest.mark.unit
def test_sample_config_structure(sample_config):
    """Test that sample config has required fields"""
    required_fields = [
        'projectName', 'dataDir', 'outDir', 'sampleFile',
        'fq_fwd', 'fq_rvr'
    ]
    for field in required_fields:
        assert field in sample_config, f"Missing required field: {field}"


@pytest.mark.unit
def test_config_files_exist(repo_root):
    """Test that example config files exist"""
    config_dir = repo_root / "configs"
    assert config_dir.exists()

    # Check for test configs
    test_configs = [
        "test_variant_calling_config.yaml",
        "test_assembly_config.yaml",
        "basic_config.yaml"
    ]

    for config_file in test_configs:
        config_path = config_dir / config_file
        assert config_path.exists(), f"Missing config file: {config_file}"


@pytest.mark.unit
def test_config_yaml_valid(repo_root):
    """Test that config files are valid YAML"""
    config_dir = repo_root / "configs"

    for config_file in config_dir.glob("*.yaml"):
        try:
            with open(config_file) as f:
                yaml.safe_load(f)
        except yaml.YAMLError as e:
            pytest.fail(f"Invalid YAML in {config_file}: {e}")
