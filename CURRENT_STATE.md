# NCCR Genomics Pipeline - Current State Documentation

**Date:** 2026-01-26
**Branch:** refactor-snakemake9-python312
**Purpose:** Comprehensive inventory before refactoring

---

## Current Version & Dependencies

### Python & Snakemake Versions
- **Current Python:** 3.8.* (as per environment.yaml)
- **Current Snakemake:** >=7.1.* (as per environment.yaml)
- **README Claims:** Snakemake 5.22 ⚠️ **INCONSISTENT**
- **Target Python:** 3.12.*
- **Target Snakemake:** 9.x

### Dependencies
- click (no version specified)
- pyaml (no version specified)

---

## CLI Commands Available

### Core Workflow Commands

1. **`samples`** - Generate samplesheet from FastQ directory
   - Status: ✅ Active utility
   - Uses: fastq_dir_to_samplesheet script

2. **`clean`** - Run preprocessing workflow
   - Status: ✅ Active
   - Snakemake target: `preprocess`

3. **`assemble`** - Run assembly pipeline
   - Status: ✅ Active
   - Snakemake target: `assemble`

4. **`align`** - Run alignment workflow
   - Status: ✅ Active
   - Snakemake target: `align`

5. **`call`** - Run breseq variant calling
   - Status: ✅ Active
   - Snakemake target: `breseq`

6. **`funcall`** - Run fungal variant calling
   - Status: ✅ Active
   - Snakemake target: `varcall`

7. **`annotate`** - Run annotation pipeline
   - Status: ✅ Active
   - Snakemake target: `annotate`

8. **`gapseq`** - Run Gapseq pipeline
   - Status: 🧪 Experimental (comment says "New and exploratory")
   - Snakemake target: `run_gapseq`

9. **`isolate`** - Main isolate genomics pipeline
   - Status: ✅ **PRIMARY WORKFLOW**
   - Methods: `call_variants`, `assemble`, `assemble_only`
   - Most flexible command with method selection

10. **`unlock`** - Unlock Snakemake working directory
    - Status: ✅ Utility command

11. **`metagenome`** - Metagenomic assembly comparison
    - Status: 🚧 Under construction (comment in code)
    - Uses: Snakefile_metagenome
    - Method: `metaflye`

12. **`phage`** - Phage detection with geNomad
    - Status: ✅ Active (uses Snakefile_test)
    - Snakemake target: `find_phage`

### ❌ Missing RNAseq Command
- README documents RNAseq workflows (`star`, `kallisto`)
- **NO CLI command found in main.py**
- Config examples exist: `configs/rnaseq_config.yaml`
- Rules exist: mentioned in README but not linked

---

## Snakemake Rule Files

### Active Rules (12 files)
1. **preprocess.smk** - Read QC and preprocessing
2. **preprocess_old.smk** - ⚠️ Old version (to be removed)
3. **assemble.smk** - Genome assembly (SPAdes/Unicycler)
4. **annotate.smk** - Gene annotation (Prokka, eggNOG)
5. **call_variants.smk** - Variant calling (BWA, breseq, bcftools)
6. **alignment.smk** - Read alignment
7. **compare_genomes.smk** - Genome comparisons
8. **count.smk** - Read counting
9. **profile.smk** - Profiling tools
10. **typing.smk** - Strain typing
11. **quast.smk** - Assembly QC
12. **hybrid_assembly.smk** - Hybrid assembly (PacBio/Illumina)

---

## Conda Environments (16 files)

1. **qc.yaml** - Preprocessing (BBMap, FastQC)
2. **assemble.yaml** - Assembly tools (SPAdes, etc.)
3. **unicycler.yaml** - Unicycler assembler
4. **hybrid_assembly.yaml** - Hybrid assembly tools
5. **annotate.yaml** - Prokka annotation
6. **emapper.yaml** - eggNOG-mapper
7. **call_variants.yaml** - Variant calling (BWA, bcftools, breseq)
8. **align.yaml** - Alignment tools
9. **anVar.yaml** - Variant annotation (SnpEff)
10. **compare_genomes.yaml** - FastANI, MUMmer
11. **quast.yaml** - QUAST
12. **typing.yaml** - MLST, serotyping
13. **count.yaml** - Read counting
14. **profile.yaml** - Profiling tools
15. **orthoFinder.yaml** - OrthoFinder
16. **panX.yaml** - Pan-genome analysis

---

## Snakefile Variants

1. **Snakefile** - Main workflow (isolate genomics)
2. **Snakefile_metagenome** - Metagenomic workflows
3. **Snakefile_test** - Test/phage workflows

---

## Workflow Capabilities

### Preprocessing
- BBMap trimming and filtering
- Adapter removal
- PhiX removal
- Host decontamination (human/mouse)
- FastQC quality control
- Optional read merging

### Assembly
- **SPAdes** (default)
- **Unicycler** (optional, modes: conservative/normal/bold)
- **Hybrid assembly** (short + long reads)
- Assembly QC with QUAST

### Annotation
- **Prokka** - Gene calling and annotation
- **eggNOG-mapper** - Functional annotation
- **geNomad** - Phage/viral detection

### Variant Calling
- **breseq** - Mutation detection (polymorphic mode)
- **BWA + bcftools** - Traditional variant calling
- Duplicate removal with samtools
- Variant filtering
- SnpEff annotation

### Comparative Genomics
- **FastANI** - Average Nucleotide Identity
- **MUMmer/nucmer** - Genome alignment
- **PhyloPhlAn** - Phylogenetics
- **OrthoFinder** - Orthology analysis
- **PanX** - Pan-genome analysis

### Typing & Characterization
- **MLST** - Multi-locus sequence typing
- **Serotyping**
- **Gapseq** - Metabolic pathway prediction (experimental)

### Metagenomic Tools
- **mOTUs** - Taxonomic profiling
- **fetchMG** - Marker gene extraction
- **Mash** - Distance estimation

---

## Configuration Files

### Test Configs
- `configs/test_variant_calling_config.yaml` - Default test config
- `configs/test_assembly_config.yaml`
- `configs/basic_config.yaml`

### Project Configs (examples from real runs)
- Multiple dated configs (e.g., `15-02-2023-Erec-config.yaml`)
- Organism-specific configs (Btheta, Erec, oligoMM, rhizopus)

### RNAseq Configs
- `configs/rnaseq_config.yaml` - Example RNAseq config
- `configs/rnaseq_samples.txt`

---

## Test Data

Located in: `workflow/test_data/varcall_test_data/`
- Reference genome: `LL6_1.fasta.gz`
- Sample: `LL23` (paired-end reads)
- Sample file: `test_samples.txt`

---

## Python Scripts

### Workflow Scripts (in workflow/scripts/)
1. **configure_project.py** - Project setup
2. **fastq_dir_to_samplesheet.py** - Sample sheet generation
3. **contig_filter.py** - Contig filtering
4. **genlen.py** - Genome length calculation
5. **get_vars.py** - Variant extraction utilities
6. **parse_orthologs.py** - Orthology parsing
7. **quality_control_report.py** - QC reporting
8. **saturation_curve.py** - Sequencing saturation
9. **sequence.py** - Sequence utilities
10. **skew.py** - GC skew analysis
11. **snpEff_db.py** - SnpEff database management
12. **strain_var.py** - Strain variation analysis
13. **variant_base_quality.py** - Variant quality analysis

### Rule Scripts (in workflow/rules/scripts/)
1. **hifi_or_clr.py** - PacBio read type detection
2. **scratch_pad.py** - ⚠️ Development artifact (to be removed)

---

## Known Issues & Technical Debt

### Critical Issues
1. **Version inconsistency** - README vs environment.yaml (Snakemake version)
2. **Security issue** - `yaml.load()` should be `yaml.safe_load()` (main.py:42)
3. **Missing RNAseq CLI** - Documented but no command available

### Code Quality Issues
1. **Commented code** - Extensive commented blocks in main.py (lines 262-323) and Snakefile
2. **Scratch files** - Multiple `scratch_pad.py` files
3. **Old files** - `preprocess_old.smk` still present
4. **Hardcoded paths** - Some configs have hardcoded cluster paths
5. **Inconsistent echo messages** - Many commands say "Running Assembly Pipeline" incorrectly
6. **Copy-paste errors** - Multiple functions have identical implementations

### Cluster Execution
- Currently: Hardcoded SLURM submission in main.py
- Old SGE code commented out
- Hardcoded partition: "institute"
- No support for different cluster types via config
- Should use: Snakemake profiles

### Documentation
- README marked "Under Development"
- Incomplete feature documentation
- No migration guides
- Hardcoded institutional paths in examples

---

## Active vs Experimental Features

### ✅ Production-Ready
- Isolate genome assembly (SPAdes)
- Variant calling (breseq, bcftools)
- Preprocessing (BBMap)
- Basic annotation (Prokka)
- Phage detection (geNomad)

### 🧪 Experimental
- Gapseq (marked as "New and exploratory")
- Hybrid assembly (less tested)
- Pan-genome analysis
- Some typing features

### 🚧 Under Construction
- Metagenome workflows
- RNAseq CLI integration
- Many features listed in README without clear status

### ❓ Unclear Status
- PhyloPhlAn integration
- ARIBA antibiotic resistance
- Some comparative genomics features

---

## Dependencies on External Resources

### Required Databases
- Adapter sequences: `data/adapters/adapters.fa`
- PhiX genome: `data/adapters/phix174_ill.ref.fa.gz`

### Optional Databases (for specific features)
- BBMap human/mouse references (hardcoded institutional paths)
- SnpEff databases (created on-the-fly)
- eggNOG databases
- geNomad databases

---

## Cluster Integration

### Current Setup
- **Cluster Type:** SLURM
- **Default Partition:** institute
- **Job Management:** Direct submission via sbatch
- **Resource Specification:** Via `params:` in rules (mem, time, threads)

### Issues
- No support for different clusters without code modification
- Hardcoded partition name
- No resource profiles
- Old SGE code still present (commented)

---

## Backwards Compatibility Concerns

### Users Currently On
- Likely Snakemake 7.x (based on environment.yaml)
- Python 3.8
- SLURM cluster at ETH Zurich

### Breaking Changes in Migration
1. **Snakemake 9:** API changes, conda environment handling
2. **Python 3.12:** Potential deprecated syntax
3. **Cluster execution:** Moving to profiles
4. **File structure:** Reorganization may break paths

---

## Next Steps (Phase 1 Remaining)

1. ✅ Branch created: `refactor-snakemake9-python312`
2. ✅ Current state documented
3. 🔜 Set up testing infrastructure
4. 🔜 Resolve version inconsistencies
5. 🔜 Document Snakemake 7→9 API changes

---

## Questions for User

1. **Which workflows are actively used in production?**
   - Isolate variant calling?
   - Assembly?
   - RNAseq?
   - Metagenomics?

2. **Can we drop experimental features?**
   - Gapseq?
   - Pan-genome tools?
   - PhyloPhlAn?

3. **What Snakemake version is currently working?**
   - README says 5.22
   - Environment.yaml says >=7.1
   - What's actually deployed?

4. **Cluster requirements:**
   - Stay with SLURM only?
   - Need multi-cluster support?
   - Partition name configurable or always "institute"?

5. **RNAseq status:**
   - Should we restore the RNAseq CLI command?
   - Is it actively used?
   - Or can it be deprecated?
