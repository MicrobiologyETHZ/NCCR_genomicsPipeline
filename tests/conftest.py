"""
Shared pytest fixtures for NCCR Genomics Pipeline tests
"""
import pytest
from pathlib import Path
import tempfile
import shutil


@pytest.fixture
def repo_root():
    """Return the repository root directory"""
    return Path(__file__).parent.parent


@pytest.fixture
def test_data_dir(repo_root):
    """Return the test data directory"""
    return repo_root / "workflow" / "test_data" / "varcall_test_data"


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs"""
    temp = tempfile.mkdtemp(prefix="nccrpipe_test_")
    yield Path(temp)
    # Cleanup after test
    shutil.rmtree(temp, ignore_errors=True)


@pytest.fixture
def sample_config(test_data_dir, temp_dir):
    """Create a minimal test configuration"""
    return {
        'projectName': 'Test_Project',
        'dataDir': str(test_data_dir / 'raw'),
        'outDir': str(temp_dir / 'output'),
        'sampleFile': str(test_data_dir / 'test_samples.txt'),
        'fq_fwd': '_R1.fq.gz',
        'fq_rvr': '_R2.fq.gz',
        'qc': 'yes',
        'mink': 11,
        'trimq': 14,
        'mapq': 20,
        'minlen': 45,
        'merged': False,
        'fastqc': 'no'
    }


@pytest.fixture
def sample_list(test_data_dir):
    """Read and return sample list from test data"""
    sample_file = test_data_dir / 'test_samples.txt'
    if sample_file.exists():
        with open(sample_file) as f:
            return [line.strip() for line in f if line.strip()]
    return []
