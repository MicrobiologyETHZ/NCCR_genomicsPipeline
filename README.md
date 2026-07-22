# NCCR Genomics Pipeline

Pipeline for analysis of isolate genomes.

### Steps Currently Included

- Preprocessing (BBMap: adapter trimming, PhiX removal, quality filtering)
- Isolate genome assembly (SPAdes or Unicycler)
- Gene calling and functional annotation (Prokka, eggNOG-mapper)
- Phage prediction and annotation (geNomad, Cenote-Taker 3, pharokka, phold, phynteny, CheckV)
- Variant calling: **Breseq** (mutation detection) and **bcftools** (SNP calling)
- Strain-level microdiversity: **InStrain** (profile + cross-sample compare)


## Installation

Requires `conda` and Python >= 3.12.

```bash
git clone https://github.com/MicrobiologyETHZ/NCCR_genomicsPipeline
cd NCCR_genomicsPipeline
conda env create -f environment.yaml
conda activate nccrPipe
pip install -e .
```

Conda environments for individual pipeline steps are in `workflow/envs/` and are NOT created automatically. Create them first, naming each env after its YAML file. For example, the InStrain step uses the `instrain` env:

```bash
conda env create -n instrain -f workflow/envs/instrain.yaml
```


## Running Breseq

### 1. Prepare input files

**Raw reads** — paired-end FASTQ files, one directory per sample:

```
dataDir/
├── Sample1/
│   ├── Sample1_R1.fq.gz
│   └── Sample1_R2.fq.gz
└── Sample2/
    ├── Sample2_R1.fq.gz
    └── Sample2_R2.fq.gz
```

**Sample sheet** — CSV file with columns `sample`, `unit`, `fastq_1`, `fastq_2`:

```
sample,unit,fastq_1,fastq_2
Sample1,1,/path/to/dataDir/Sample1/Sample1_R1.fq.gz,/path/to/dataDir/Sample1/Sample1_R2.fq.gz
Sample2,1,/path/to/dataDir/Sample2/Sample2_R1.fq.gz,/path/to/dataDir/Sample2/Sample2_R2.fq.gz
```

Use the `samples` command to generate this file automatically from your FastQ directory (see [Generating the sample sheet](#generating-the-sample-sheet) below).

**Reference genome** — GenBank (`.gbk`) file is required for breseq. If you only have a FASTA, you can also supply a separate GFF annotation (see config options below).


### 2. Generate the sample sheet

The `samples` command scans a directory of FastQ files and writes the CSV sample sheet for you.

**Option A — from a config file** (recommended: reads `dataDir`, `samples`, `fq_fwd`, `fq_rvr`, and name-sanitising settings directly from config):

```bash
nccrPipe samples --configfile /path/to/your_config.yaml
```

**Option B — from command-line flags:**

```bash
nccrPipe samples -i /path/to/dataDir -o /path/to/samples.csv -r1 _R1.fq.gz -r2 _R2.fq.gz
```

#### Sample name sanitisation

By default the sample name is derived from the FastQ filename by stripping the read extension. If your filenames contain extra metadata (e.g. `Sample1_LIBXXX_R1.fq.gz`) you can trim them down using the sanitise options:

| Option | Config key | Default | Description |
|--------|-----------|---------|-------------|
| `-sn` / `--sanitise_name` | `sanitise_name: true` | off | Enable name sanitisation |
| `-sd` / `--sanitise_name_delimiter` | `name_delimiter: _` | `_` | Split filename on this character |
| `-si` / `--sanitise_name_index` | `name_index: 1` | `1` | Keep only the first N parts (1-based) |

**Example:** filename `Sample1_LIBXXX_R1.fq.gz` with `name_delimiter: _` and `name_index: 1` → sample name `Sample1`.

These settings live in the config file (see template) and are also used when running `samples --configfile`.


### 4. Create a config file

Copy `configs/breseq_config.yaml` and fill in your paths. The key sections are:

```yaml
dataDir: /path/to/raw/fastq/files
outDir: /path/to/output
samples: /path/to/samples.csv

# Adapter/PhiX references for preprocessing
adapters: /path/to/adapters.fa
phix: /path/to/phix174_ill.ref.fa.gz

# Reference genome
reference:
  refgbk: /path/to/reference.gbk

# Breseq options (see modes below)
breseq:
  polymorphic: true
  limit_fold_coverage: 200
  min_mapping_quality: 30
```


### 5. Breseq modes

All breseq options are controlled via the `breseq:` block in the config. The pipeline builds the command from these settings — no need to edit any pipeline files.

#### Clonal mode (default consensus calling)

Use when sequencing a clonal isolate. Detects fixed mutations relative to the reference.

```yaml
breseq:
  polymorphic: false
```

Runs: `breseq -j 8 -o {out_dir} -r {reference.gbk} reads_R1.fq.gz reads_R2.fq.gz`

#### Polymorphic mode

Use when sequencing a mixed population or when you want to detect variants at any frequency.

```yaml
breseq:
  polymorphic: true
```

Runs: `breseq -p -j 8 -o {out_dir} -r {reference.gbk} reads_R1.fq.gz reads_R2.fq.gz`

#### Polymorphic mode with coverage limit and mapping quality filter

```yaml
breseq:
  polymorphic: true
  limit_fold_coverage: 200     # -l: downsample reads above this fold-coverage
  min_mapping_quality: 30      # -m: minimum mapping quality to use a read
```

Runs: `breseq -p -j 8 -l 200 -m 30 -o {out_dir} -r {reference.gbk} reads_R1.fq.gz reads_R2.fq.gz`

#### Using a FASTA + GFF reference (instead of GenBank)

If you don't have a `.gbk` file, provide a FASTA genome and a separate GFF annotation:

```yaml
reference:
  refgbk: /path/to/reference.fasta   # FASTA genome (used as primary -r)
  refgff: /path/to/reference.gff     # GFF annotation (added as second -r)

breseq:
  polymorphic: true
  read_length: 200
  min_coverage: 30
```

Runs: `breseq -p -j 8 -l 200 -m 30 -o {out_dir} -r reference.fasta -r reference.gff reads_R1.fq.gz reads_R2.fq.gz`

> **Note:** Omitting `limit_fold_coverage` or `min_mapping_quality` (or setting them to `null`) will skip the corresponding flag and use breseq's built-in defaults.


### 6. Run the pipeline

Always do a dry run first to check that inputs resolve correctly:

```bash
nccrPipe call --config /path/to/your_config.yaml --dry
```

Then run:

```bash
# On a SLURM cluster
nccrPipe call --config /path/to/your_config.yaml

# On a local machine
nccrPipe call --config /path/to/your_config.yaml --local
```


### 7. Outputs

Breseq output for each sample is written to:

```
outDir/
└── breseq/
    ├── Sample1/
    │   ├── output/index.html      # Main HTML report
    │   ├── output/summary.html    # Mutation summary
    │   └── data/output.vcf        # Variants in VCF format
    └── Sample2/
        └── ...
```

Open `output/index.html` in a browser to view results.


### Troubleshooting

- Log files for each step are in `outDir/logs/`
- Run with `--dry` first to catch config errors before submitting jobs
- Breseq requires a reference with annotation (GenBank preferred). FASTA-only references will not produce gene-level output.


## Running InStrain

[InStrain](https://instrain.readthedocs.io/) measures strain-level microdiversity (nucleotide diversity, SNVs, popANI) from reads mapped to a reference, and compares strains across samples. It reuses the same alignment step as variant calling — reads are mapped to your reference FASTA with BWA, then `inStrain profile` runs per sample and (optionally) `inStrain compare` runs across samples.

### 1. Prepare input files

Same preprocessing inputs as Breseq (sample sheet + raw reads). InStrain additionally needs a **reference FASTA** (not GenBank): the file `<refDir>/<refName>.fasta`.

### 2. Create a config file

Copy `configs/instrain_config.yaml` and fill in your paths. The key sections:

```yaml
# Reference(s): refDir holds each <name>.fasta; <name> maps to the samples
# mapped against it. References with >=2 samples are eligible for `compare`.
reference:
  refDir: /path/to/refs
  myref:
    - Sample1
    - Sample2

instrain:
  min_cov: 5          # --min_cov for profiling
  run_compare: true   # also run `inStrain compare` across samples
  # genes:            # OPTIONAL: prodigal-style gene calls for gene-level stats
  #   myref: /path/to/myref.genes.fna
```

To get gene-level microdiversity (per-gene coverage, dN/dS), supply a prodigal-style
gene-calls FASTA per reference under `instrain.genes` (keyed by reference name); it is
passed to InStrain via `-g`. Omit it for scaffold-level profiling only.

### 3. Run the pipeline

```bash
# Per-sample profiling (dry run first)
nccrPipe instrain --config /path/to/your_config.yaml --dry
nccrPipe instrain --config /path/to/your_config.yaml          # cluster
nccrPipe instrain --config /path/to/your_config.yaml --local  # local machine

# Also run cross-sample comparison (popANI / strain sharing)
nccrPipe instrain --config /path/to/your_config.yaml --compare
```

Runs (per sample): `inStrain profile {bam} {ref}.fasta -o {out} -p {threads} [-g genes.fna] --min_cov {min_cov}`
Runs (per reference, ≥2 samples): `inStrain compare -i {profile1} {profile2} ... -o {out} -p {threads}`

### 4. Outputs

```
outDir/
├── instrain/
│   ├── profiles/
│   │   ├── Sample1_to_myref/      # per-sample IS profile (genome_info.tsv, SNVs.tsv, ...)
│   │   └── Sample2_to_myref/
│   └── compare/
│       └── myref/                 # cross-sample comparison (comparisonsTable.tsv, ...)
└── logs/instrain/                 # per-step logs
```


## Running Genome Assembly

### 1. Prepare input files

**Raw reads** — paired-end FASTQ files (one pair per sample, anywhere on disk):

```
/path/to/raw/
├── Sample1_R1.fq.gz
├── Sample1_R2.fq.gz
├── Sample2_R1.fq.gz
└── Sample2_R2.fq.gz
```

**Sample sheet** — CSV with columns `sample`, `unit`, `fastq_1`, `fastq_2`. Generate automatically:

```bash
nccrPipe samples -c configs/assembly_config.yaml
```

### 2. Create a config file

Copy `configs/assembly_config.yaml` and fill in your paths:

```yaml
data_dir: /path/to/raw/fastq
output_dir: /path/to/output
samples: /path/to/samples.csv

adapters: data/adapters/adapters.fa     # bundled in the repo
phix: data/adapters/phix174_ill.ref.fa.gz

assembler: spades   # or unicycler
```

### 3. Run assembly (SPAdes + Prokka)

```bash
# Dry run first
nccrPipe assemble -c configs/assembly_config.yaml --dry

# Local machine
nccrPipe assemble -c configs/assembly_config.yaml --local

# SLURM cluster
nccrPipe assemble -c configs/assembly_config.yaml
```

### 4. Run functional annotation (eggNOG-mapper)

Requires the eggNOG database.

Add to your config:

```yaml
database: eggnog
eggnog_db: /path/to/eggnog-data
```

Then run:

```bash
nccrPipe annotate -c configs/assembly_config.yaml --local
```

### 5. Outputs

```
output_dir/
├── clean_reads/{sample}/          # QC-filtered reads
├── assembly/{sample}/
│   ├── {sample}.scaffolds.min200.fasta   # assembled scaffolds
│   ├── {sample}.gff                      # Prokka annotation
│   └── {sample}.faa                      # protein sequences
└── logs/
```


## Running Phage Prediction and Annotation

Predicts phages/proviruses with **geNomad** and **Cenote-Taker 3**, then annotates
each caller's viral contigs through **pharokka → phold → phynteny_transformer**.
**CheckV** assesses completeness and contamination. The two callers are annotated
separately so their results can be compared directly.

Unlike the rest of the pipeline, this step is not driven by `samples.csv` — it is
usually pointed at assemblies produced elsewhere (public genomes, collaborator
data), and works fine on read-only input directories.

### 1. Install the databases (once)

Six databases, ~20 GB total. Create the tool environments, then fetch them:

```bash
for e in genomad cenotetaker checkv pharokka phold phynteny; do
    conda env create -n "$e" -f workflow/envs/"$e".yaml
done

bash workflow/scripts/setup_phage_dbs.sh /path/to/phage_dbs
```

The script prints the resulting paths for your config. Note CheckV's directory
carries a version suffix (`checkv-db-v1.5`) that changes between releases.

### 2. Create a config file

Start from `configs/phage_config.yaml`. Assemblies can be supplied three ways,
and the forms can be combined:

```yaml
outDir: ../output

# 1. A directory of assemblies (the usual case). Gzipped files are handled
#    automatically; *.fa / *.fasta / *.fna / *.fas are matched by default.
assembly_dir: /path/to/assemblies
# pattern: "*.fna"

# 2. Specific files — a list, or a name -> path map when you need to control
#    output names (required if two assemblies share a filename).
# assemblies:
#   strainA: /path/to/runA/scaffolds.fasta
#   strainB: /path/to/runB/scaffolds.fasta

# 3. Assemblies built by this pipeline.
# samples: samples.csv
# assembler: spades

callers: [genomad, cenotetaker]   # either or both
annotate: true                    # false = prediction only

databases:
  genomad: /path/to/phage_dbs/genomad_db
  cenotetaker: /path/to/phage_dbs/ct3_dbs
  checkv: /path/to/phage_dbs/checkv-db-v1.5
  pharokka: /path/to/phage_dbs/pharokka_db
  phold: /path/to/phage_dbs/phold_db
  phynteny: /path/to/phage_dbs/phynteny_models
```

Assemblies are named by filename stem (`GCA_000001_genomic.fna.gz` →
`GCA_000001_genomic`). If two files share a stem — several runs each producing
`scaffolds.fasta` — the parent directory is prepended. If that still collides,
the workflow stops and asks for explicit names rather than overwriting results.

Only the databases the enabled steps need are checked, so a prediction-only run
(`annotate: false`) does not require the annotation databases.

By default the rules use the pinned `workflow/envs/*.yaml` files, so `--use-conda`
builds exactly the versions this workflow targets. To reuse environments you have
already installed elsewhere, point at their parent directory instead:

```yaml
conda_env_dir: /nfs/.../conda_envs   # expects genomad/, checkv/, phold/, ... inside
```

### 3. Run the pipeline

```bash
# Dry run first
nccrPipe phage -c /path/to/phage_config.yaml --dry

# Prediction + annotation + summary tables
nccrPipe phage -c /path/to/phage_config.yaml

# Prediction only
nccrPipe phage -c /path/to/phage_config.yaml --predict-only

# Local machine
nccrPipe phage -c /path/to/phage_config.yaml --local
```

phold's structure prediction runs on CPU by default and is the slowest step; set
`phold: {cpu: false}` only if submitting to a GPU partition.

### 4. Outputs

```
output_dir/phage/
├── {assembly}/
│   ├── genomad/                                   # geNomad end-to-end output
│   ├── cenotetaker/                               # Cenote-Taker 3 output
│   ├── viral/{assembly}.{caller}.fna              # predicted viral contigs
│   ├── checkv/{caller}/quality_summary.tsv        # completeness/contamination
│   └── annotate/{caller}/
│       ├── pharokka/pharokka.gbk
│       ├── phold/phold.gbk
│       └── phynteny/phynteny_transformer.gbk      # final annotation
└── summary/
    ├── phage_predictions.tsv                      # one row per viral contig
    └── phage_annotations.tsv                      # one row per assembly x caller
```

Assemblies with no phage are normal and do not fail the run: the viral FASTA is
empty, the annotation steps are skipped, and the assembly appears in the summary
tables with zero contigs.

### Version note

Bioconda has no geNomad 1.9.1 (releases go 1.9.0 → 1.10.0), so `envs/genomad.yaml`
pins 1.12.0. If you are reproducing a collaborator's results, confirm which
version they actually ran.


## Running Other Pipeline Steps

```bash
# Preprocessing only
nccrPipe isolate -c /path/to/config.yaml -m preprocess --local

# BCFtools variant calling (SNP calling against reference FASTA)
nccrPipe isolate -c /path/to/config.yaml -m call_variants
```
