# NCCR Genomics Pipeline - Test Suite

## Overview

This directory contains the test suite for the NCCR Genomics Pipeline refactoring project.

## Test Structure

```
tests/
├── unit/              # Unit tests for individual functions
├── integration/       # Integration tests for workflows
├── fixtures/          # Shared test data and fixtures
└── conftest.py        # Shared pytest fixtures
```

## Running Tests

### Install Test Dependencies

```bash
pip install pytest pytest-cov
```

### Run All Tests

```bash
pytest
```

### Run Specific Test Categories

```bash
# Run only unit tests
pytest -m unit

# Run only integration tests
pytest -m integration

# Run tests excluding slow tests
pytest -m "not slow"

# Run tests excluding cluster-dependent tests
pytest -m "not cluster"
```

### Run with Coverage

```bash
pytest --cov=workflow --cov-report=html
```

## Test Markers

Tests are marked with the following categories:

- `unit`: Unit tests for individual functions
- `integration`: Integration tests for complete workflows
- `slow`: Tests that take a long time to run
- `cluster`: Tests that require cluster access
- `data`: Tests that require test data

## Writing Tests

### Unit Test Example

```python
import pytest

@pytest.mark.unit
def test_my_function():
    result = my_function(input_data)
    assert result == expected_output
```

### Integration Test Example

```python
import pytest

@pytest.mark.integration
@pytest.mark.data
def test_workflow(temp_dir, sample_config):
    # Test complete workflow
    pass
```

## Continuous Integration

Tests run automatically on GitHub Actions for:
- All pushes to master/main/refactor branches
- All pull requests

See `.github/workflows/tests.yml` for CI configuration.

## Test Data

Test data is located in `workflow/test_data/varcall_test_data/` and includes:
- Sample: LL23 (paired-end reads)
- Reference genome: LL6_1.fasta.gz
- Sample list: test_samples.txt

## TODO

- [ ] Add more unit tests for Python scripts
- [ ] Add integration tests for each workflow type
- [ ] Add tests for CLI commands
- [ ] Add tests for config validation
- [ ] Add performance benchmarks
