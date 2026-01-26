# Refactoring Decisions - NCCR Genomics Pipeline

**Date:** 2026-01-26
**Stakeholder:** Anna Sintsova
**Purpose:** Document decisions to guide refactoring scope

---

## Current Working Environment

- **Snakemake Version:** 7.32.4 (CONFIRMED)
- **Python Version:** 3.8
- **Target:** Snakemake 9.x, Python 3.12

---

## Features to KEEP (Active Production Use)

### Core Workflows
1. ✅ **Isolate variant calling**
   - breseq (polymorphic mode)
   - bcftools variant calling
   - BWA alignment
   - Priority: HIGH

2. ✅ **Eukaryotic variant calling**
   - Generalized from fungal-specific (funcall)
   - BWA + bcftools for eukaryotes
   - Priority: MEDIUM

3. ✅ **Genome assembly (ALL methods needed)**
   - SPAdes (isolate, primary)
   - Unicycler (isolate, secondary)
   - **Hybrid assembly** (PacBio + Illumina) - KEEP
   - **Long read metagenome assembly** - KEEP
   - Priority: HIGH

4. ✅ **Preprocessing**
   - BBMap quality control
   - Adapter removal
   - PhiX filtering
   - **Host decontamination (human/mouse)** - Make configurable
   - Priority: HIGH

5. ✅ **Annotation**
   - Prokka (gene calling)
   - eggNOG-mapper (functional annotation)
   - **Gapseq** (EXPAND - metabolic pathway prediction)
   - Priority: HIGH

5. ✅ **Phage detection**
   - geNomad
   - Priority: MEDIUM

### Comparative/Analysis Tools
6. ✅ **mOTUs**
   - Taxonomic profiling
   - Priority: MEDIUM

7. ✅ **OrthoFinder**
   - Orthology analysis
   - Priority: MEDIUM

### Utilities
8. ✅ **FastQC** - Quality control reporting
9. ✅ **QUAST** - Assembly QC
10. ✅ **FastANI** - Average Nucleotide Identity (if actively used)
11. ✅ **MUMmer/nucmer** - Genome comparison (if actively used)

---

## Features to REMOVE (Not Used)

### Complete Removal
1. ❌ **RNAseq workflows**
   - STAR pipeline
   - kallisto pipeline
   - featureCounts
   - All RNAseq rules, configs, documentation
   - Reason: Separate pipeline exists
   - Files to remove:
     - `configs/rnaseq_config.yaml`
     - `configs/rnaseq_samples.txt`
     - `workflow/rules/count.smk`
     - `workflow/rules/alignment.smk` (if RNAseq specific)
     - `workflow/envs/count.yaml`
     - RNAseq documentation in README

2. ⚠️ **Metagenomics workflows (PARTIAL REMOVAL)**
   - ❌ Remove: Short read metagenomics workflows
   - ❌ Remove: Snakefile_metagenome (if not for long read assembly)
   - ✅ KEEP: Long read metagenome assembly (PacBio + Illumina)
   - ✅ KEEP: Hybrid assembly rules
   - Reason: Long read assembly needed, other metagenomics not used

3. ❌ **Pan-genome analysis**
   - PanX
   - Reason: Not needed
   - Files: `workflow/envs/panX.yaml`, rules in Snakefile

4. ❌ **Antibiotic resistance**
   - ARIBA
   - Reason: Not needed

5. ❌ **Strain typing**
   - MLST
   - Serotyping
   - Files: `workflow/rules/typing.smk`, `workflow/envs/typing.yaml`

6. ❌ **Phylogenetics**
   - PhyloPhlAn
   - Reason: Not needed
   - Remove from documentation and rules

### CLI Commands to Remove/Refactor
- `clean` ❌ Remove (redundant with isolate -m preprocess)
- `align` ❌ Remove (redundant with isolate workflow)
- `call` ❌ Remove (redundant with isolate -m call_variants)
- `funcall` ✏️ **Rename to `eukaryote`** (generalize fungal → eukaryotic)
- `annotate` ❌ Remove (redundant with isolate workflow)
- `gapseq` ❌ Remove command (integrate into isolate workflow)
- `assemble` ❌ Remove (redundant with isolate -m assemble)
- `unlock` ✅ Keep (utility)
- `metagenome` ❌ Remove (long read assembly via isolate with method)
- `phage` ✅ Keep

**Streamline to:**
- `samples` - Generate sample sheets (utility)
- `isolate` - Main workflow with methods:
  - preprocess
  - assemble (SPAdes/Unicycler)
  - hybrid_assemble (PacBio + Illumina)
  - call_variants (bacterial)
  - annotate (Prokka/eggNOG/Gapseq)
- `eukaryote` - Eukaryotic variant calling
- `phage` - Phage detection
- `profile` - mOTUs profiling
- `orthology` - OrthoFinder
- `unlock` - Utility

---

## Cluster Configuration Requirements

### Multi-Cluster Support
- ✅ **SLURM** (primary)
- ✅ **SGE** (keep support)
- ✅ Other queueing systems via profiles

### Configurable Parameters
- ✅ **Partition name** - User configurable (institute, sunagawa, etc.)
- ✅ **Queue name** - For different systems
- ✅ **Resource defaults** - Configurable per rule

### Implementation
- Use Snakemake profiles for flexibility
- Create templates for common clusters:
  - `profiles/slurm-institute/`
  - `profiles/slurm-sunagawa/`
  - `profiles/sge/`
  - `profiles/local/`

---

## File Structure Reorganization

### Current Issues
- setup.py references 'workflow' as package (incorrect)
- Configs directory has many test/project configs
- Multiple Snakefile variants

### Proposed Structure
```
NCCR_genomicsPipeline/
├── nccrpipe/              # Python package
│   ├── __init__.py
│   ├── cli.py             # CLI commands
│   ├── utils.py           # Utilities
│   └── scripts/           # Python scripts
├── workflow/              # Snakemake workflows
│   ├── Snakefile          # Main workflow
│   ├── rules/             # Rule modules
│   │   ├── preprocess.smk
│   │   ├── assemble.smk
│   │   ├── call_variants.smk
│   │   ├── annotate.smk
│   │   ├── compare_genomes.smk
│   │   └── hybrid_assembly.smk
│   ├── envs/              # Conda environments
│   ├── scripts/           # Shell/helper scripts
│   └── profiles/          # Cluster profiles
│       ├── slurm/
│       ├── sge/
│       └── local/
├── config/                # Example configs only
│   ├── example_variant_calling.yaml
│   ├── example_assembly.yaml
│   └── README.md
├── test_data/             # Move from workflow/
├── tests/                 # Test suite
├── docs/                  # Documentation
├── environment.yaml       # Main environment
├── pyproject.toml         # Modern packaging
└── README.md
```

---

## Snakefile Consolidation

### Current
- `Snakefile` - isolate workflows
- `Snakefile_metagenome` - Check if needed for long read assembly
- `Snakefile_test` - phage workflows

### Proposed
- ✅ **Single consolidated Snakefile** with modular rules
- Different target rules for different workflows
- Use `include:` statements for modularity
- All assembly methods (short, long, hybrid) in one place

---

## Priority Order for Refactoring

### Phase 2: Cleanup (NEXT)
1. **HIGH:** Remove RNAseq code completely
2. **HIGH:** Remove metagenomics code
3. **HIGH:** Remove PanX, ARIBA, typing, PhyloPhlAn
4. **MEDIUM:** Remove scratch_pad.py files
5. **MEDIUM:** Remove commented code blocks
6. **LOW:** Clean up old config files

### Phase 3: Python 3.12 Compatibility
1. **HIGH:** Test all kept scripts with Python 3.12
2. **HIGH:** Fix yaml.load() → yaml.safe_load()
3. **MEDIUM:** Update dependencies

### Phase 4: Snakemake 9 Migration
1. **CRITICAL:** Replace touch() in all kept rules
2. **CRITICAL:** Migrate params to resources
3. **HIGH:** Create cluster profiles
4. **HIGH:** Update CLI for profile-based execution

### Phase 5: Enhancements
1. **HIGH:** Integrate Gapseq into annotation workflow
2. **MEDIUM:** Streamline CLI commands
3. **MEDIUM:** Add config validation
4. **LOW:** Improve documentation

---

## Testing Strategy

### Must Test (Priority)
1. ✅ Variant calling workflow (breseq + bcftools)
2. ✅ Assembly workflow (SPAdes)
3. ✅ Preprocessing only
4. ✅ Annotation (Prokka + eggNOG)
5. ✅ Complete pipeline (preprocess → assemble → annotate → variants)

### Should Test
- Unicycler assembly
- Hybrid assembly
- Phage detection
- mOTUs profiling
- OrthoFinder

### Test Environments
- Python 3.12 + Snakemake 9.x
- Local execution
- SLURM cluster (if available)

---

## Backward Compatibility Plan

### For Current Users (Snakemake 7.32.4, Python 3.8)
- **master branch:** Keep stable with current versions
- **Tag:** Create v0.9 tag before major changes
- **Support period:** 3-6 months overlap

### Migration Path
1. **v0.9** (current) - Snakemake 7.32.4, Python 3.8
2. **v1.0-rc** (refactor branch) - Testing release
3. **v1.0** (stable) - Snakemake 9.x, Python 3.12

---

## Decisions Finalized

1. ✅ **FastANI and nucmer:** KEEP as is
2. ✅ **Hybrid assembly:** KEEP - actively used
3. ✅ **Long read metagenome assembly:** KEEP - actively used
4. ✅ **Fungal variant calling (funcall):** RENAME/GENERALIZE to eukaryotic variant calling
5. ✅ **Host decontamination paths:** Make configurable in config file
6. ✅ **Multiple Snakefiles:** Consolidate into single Snakefile
7. ✅ **All assembly methods:** Keep SPAdes, Unicycler, and hybrid assembly

---

## Summary

**Keep:** 9 core workflows (all assembly methods, variant calling, annotation, profiling)
**Remove:** 5 major features (RNAseq, PanX, ARIBA, MLST typing, PhyloPhlAn)
**Streamline:** CLI from 12 commands → 7 focused commands
**Modernize:** Cluster profiles, Python 3.12, Snakemake 9, single Snakefile
**Enhance:** Configurable host decontamination, generalized eukaryotic variant calling

**Result:** Leaner, more maintainable pipeline focused on isolate/metagenome genomics with flexible assembly
