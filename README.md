# NCCR Genomics Pipeline

Pipeline for analysis of isolate genomes.

### Steps Currently Included

- Preprocessing (BBMap: adapter trimming, PhiX removal, quality filtering)
- Isolate genome assembly (SPAdes or Unicycler)
- Gene calling and functional annotation (Prokka, eggNOG-mapper, geNomad)
- Variant calling: **Breseq** (mutation detection) and **bcftools** (SNP calling)


## Installation

Requires `conda` and Python >= 3.12.

```bash
git clone https://github.com/MicrobiologyETHZ/NCCR_genomicsPipeline
cd NCCR_genomicsPipeline
conda env create -f environment.yaml
conda activate nccrPipe
pip install -e .
```

Conda environments for individual pipeline steps are in `workflow/envs/` and are created automatically on first run via `--use-conda`.


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


## Running Other Pipeline Steps

```bash
# Preprocessing only
nccrPipe isolate -c /path/to/config.yaml -m preprocess --local

# Genome assembly
nccrPipe isolate -c /path/to/config.yaml -m assemble

# BCFtools variant calling (SNP calling against reference FASTA)
nccrPipe isolate -c /path/to/config.yaml -m call_variants
```
