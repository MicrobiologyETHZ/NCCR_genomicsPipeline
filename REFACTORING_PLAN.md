# NCCR Genomics Pipeline - Refactoring Plan

**Goal:** Update pipeline to work with Snakemake 9 and Python 3.12, clean up code, and make it usable by other people.

**Status:** Planning Phase
**Last Updated:** 2026-01-26

---

## Phase 1: Assessment & Preparation (Non-Breaking)
**Goal:** Understand current usage and prepare for safe refactoring

### 1.1 Create development/refactoring branch
- [ ] Protect master branch
- [ ] Create `refactor` or `dev` branch for all refactoring work
- [ ] Users continue using master until refactoring is complete

### 1.2 Document current state
- [ ] Inventory all active workflows (isolate, rnaseq, metagenome, phage)
- [ ] Identify which features are actively used vs. experimental
- [ ] Document current dependencies and versions
- [ ] List breaking changes in Snakemake 7→9 and Python 3.8→3.12

### 1.3 Set up testing infrastructure
- [ ] Create `tests/` directory
- [ ] Add pytest framework
- [ ] Create integration tests using existing test data
- [ ] Establish CI/CD pipeline (GitHub Actions)

### 1.4 Version inconsistency audit
- [ ] README says Snakemake 5.22, environment.yaml says >=7.1.*
- [ ] Determine actual working version
- [ ] Document Snakemake API changes needed for v9

---

## Phase 2: Code Cleanup (Mostly Non-Breaking)
**Goal:** Remove technical debt without changing functionality

### 2.1 Remove development artifacts
- [ ] Delete `scratch_pad.py` files
- [ ] Remove commented-out code blocks (extensive in main.py, Snakefile)
- [ ] Clean up unused imports
- [ ] Remove old/experimental code paths

### 2.2 Fix immediate code quality issues
- [ ] Fix deprecation warnings
- [ ] Add proper logging instead of print statements
- [ ] Standardize path handling (use pathlib consistently)
- [ ] Fix hardcoded paths in configs

### 2.3 Organize file structure
Current structure is messy. Target structure:

```
NCCR_genomicsPipeline/
├── nccrpipe/              # Python package (lowercase, no underscores)
│   ├── __init__.py
│   ├── cli.py             # CLI commands
│   ├── utils.py           # Helper functions
│   └── scripts/           # Python scripts for rules
├── workflow/              # Snakemake workflows
│   ├── Snakefile
│   ├── rules/             # Rule modules
│   ├── envs/              # Conda environments
│   └── scripts/           # Shell scripts
├── config/                # Example configs (not 'configs')
├── test_data/             # Keep test datasets
├── tests/                 # Unit & integration tests
├── docs/                  # User & developer docs
├── pyproject.toml         # Modern Python packaging
├── environment.yaml       # Main conda environment
└── README.md
```

---

## Phase 3: Python 3.12 Compatibility
**Goal:** Ensure code works with Python 3.12

### 3.1 Update Python version
- [ ] Update environment.yaml: `python==3.12.*`
- [ ] Test all Python scripts for 3.12 compatibility
- [ ] Fix any deprecated syntax/imports
- [ ] Update setup.py/pyproject.toml metadata

### 3.2 Update Python dependencies
- [ ] Update Click to latest version
- [ ] Update PyYAML to latest (yaml.load → yaml.safe_load)
- [ ] Update pandas, numpy if used
- [ ] Pin versions appropriately

---

## Phase 4: Snakemake 9 Migration
**Goal:** Update to Snakemake 9 with minimal breaking changes

### 4.1 Snakemake API updates
- [ ] Update `conda:` directive syntax (now requires full path or env name)
- [ ] Update `params:` for cluster configuration (Snakemake profiles)
- [ ] Replace deprecated `touch()` with proper output markers
- [ ] Update `configfile` handling
- [ ] Fix any `shell:` command syntax issues

### 4.2 Modernize cluster execution
- [ ] Replace hardcoded cluster commands with Snakemake profiles
- [ ] Create example SLURM profile
- [ ] Document profile usage
- [ ] Remove cluster string building from main.py

### 4.3 Update conda environments
- [ ] Test all environment.yaml files with Snakemake 9
- [ ] Update tool versions if needed
- [ ] Ensure compatibility

---

## Phase 5: Improve Package Structure
**Goal:** Make it installable and maintainable

### 5.1 Modernize packaging
- [ ] Replace setup.py with pyproject.toml (PEP 518)
- [ ] Add proper metadata (version, description, classifiers)
- [ ] Add development dependencies
- [ ] Add entry points for CLI

### 5.2 Improve CLI
- [ ] Refactor main.py into clean modules
- [ ] Add better error handling
- [ ] Add config validation
- [ ] Standardize command structure
- [ ] Add --version flag
- [ ] Improve help text

### 5.3 Configuration management
- [ ] Create config schema/validation
- [ ] Add config templates
- [ ] Better error messages for config issues
- [ ] Document all config parameters

---

## Phase 6: Documentation & Usability
**Goal:** Make it easy for others to use

### 6.1 User documentation
- [ ] Rewrite README with clear sections
- [ ] Add installation guide (conda, pip)
- [ ] Add quickstart tutorial
- [ ] Document each workflow type
- [ ] Add troubleshooting section
- [ ] Add changelog

### 6.2 Developer documentation
- [ ] Add contributing guidelines
- [ ] Document code structure
- [ ] Add docstrings to functions
- [ ] Create development setup guide
- [ ] Document testing procedures

### 6.3 Example workflows
- [ ] Provide minimal working examples
- [ ] Document expected inputs/outputs
- [ ] Add data preparation guides

---

## Phase 7: Testing & Validation
**Goal:** Ensure nothing breaks

### 7.1 Comprehensive testing
- [ ] Unit tests for Python functions
- [ ] Integration tests for workflows
- [ ] Test with actual data on each workflow type
- [ ] Performance benchmarks
- [ ] Test on different systems/clusters

### 7.2 Validation against production
- [ ] Run side-by-side with current version
- [ ] Verify outputs are identical
- [ ] Performance comparison
- [ ] User acceptance testing

---

## Phase 8: Release & Migration
**Goal:** Safe transition for users

### 8.1 Version release
- [ ] Tag v1.0.0 (breaking changes from 0.x)
- [ ] Create release notes
- [ ] Migration guide for existing users
- [ ] Deprecation warnings where appropriate

### 8.2 Documentation deployment
- [ ] Consider ReadTheDocs or GitHub Pages
- [ ] Create installation video/tutorial
- [ ] Share with user community

### 8.3 Ongoing maintenance
- [ ] Set up issue templates
- [ ] Establish release cycle
- [ ] Plan for future features

---

## Key Principles Throughout

- ✅ **Test after each phase** - Don't accumulate untested changes
- ✅ **Maintain backward compatibility when possible** - Use deprecation warnings
- ✅ **Keep master branch stable** - Users can continue working
- ✅ **Document as you go** - Don't defer documentation
- ✅ **Use branches/PRs** - Review each major change
- ✅ **Validate with real data** - Test data + production-like scenarios

---

## Timeline Estimate

- **Phases 1-2:** Foundation (1-2 weeks)
- **Phases 3-4:** Version upgrades (2-3 weeks)
- **Phases 5-6:** Improvements (2-3 weeks)
- **Phases 7-8:** Testing & release (1-2 weeks)

**Total: 6-10 weeks** depending on complexity and testing requirements

---

## Current Issues Identified

1. **Version Inconsistencies:**
   - README: Snakemake 5.22
   - environment.yaml: Snakemake >=7.1.*
   - Need to determine actual working version

2. **Code Quality:**
   - Multiple scratch_pad.py files
   - Extensive commented-out code
   - Hardcoded paths in some configs
   - Inconsistent use of Path vs string paths

3. **Package Structure:**
   - setup.py references 'workflow' package instead of proper package structure
   - No tests directory
   - Configs scattered
   - Multiple Snakefile variants (Snakefile, Snakefile_metagenome, Snakefile_test)

4. **Security Issues:**
   - yaml.load() should be yaml.safe_load()
   - Unsafe string interpolation in some shell commands

5. **Documentation:**
   - README mentions "Under Development"
   - Incomplete documentation for many features
   - No migration guide
   - No developer docs

---

## Notes

- Pipeline currently supports: isolate genomics, RNAseq, metagenomics, phage detection
- Active workflows in use: need to verify with users
- Test data available in `workflow/test_data/varcall_test_data/`
- Cluster support: SLURM (currently), previously SGE
