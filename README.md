# NCCR Genomics Pipeline

**Under Development**


Pipeline for analysis of isolate genomes. 

### Steps Currently Included:

- Preprocessing
- Isolate Genome Assembly 
- Gene calling and annotation
- Mapping back to assembly or to reference genome
- Variant calling and annotation
- ANI calculation
- Phylogenetics with PhyloPhlAn


### Installation

The software depends on `conda`, `Python >= 3.8` and `snakemake == 5.22`. 

```
git clone https://github.com/MicrobiologyETHZ/NCCR_genomicsPipeline
cd NCCR_genomicsPipeline
conda env create -f nccrPipe_environment.yaml
# mamba env create -f nccrPipe_environment.yaml
conda activate nccrPipe
pip install -e . 

```

### Conda environments

- Using pre-created environments.
- Recipes could be found in ``
- Create before running the pipeline [Write instructions]


### Running the isolate pipeline

```
Usage: nccrPipe isolate [OPTIONS]

Options:
  -c, --config TEXT  Configuration File
  -m, --method TEXT  Workflow to run, options: [call_variants, assemble,
                     assemble_only]
  --local            Run on local machine
  --no-conda         Do not use conda, under construction, do not use
  --dry              Show commands without running them
  --help             Show this message and exit.

```

### Running variant calling pipeline:

#### 1. Testing the installation.

Run ```nccrPipe isolate -m call_variants --dry```. This will do a dry run of variant calling pipeline on the test dataset. 
Run ```nccrPipe isolate -m call_variants``` This will run the variant calling pipeline on the test dataset
Use `--local` if you don't the jobs to be submitted to the queuing system. 

The first time you run, it will take a while to create the conda environments.

#### 2. Creating a config file and sample file

To run the pipeline you need to provide a `yaml` config file. Example config file for sufficient for variant calling is shown below.

```yaml
projectName: Test_Variant_Calling
# dataDir: directory with raw sequencing files, files for each sample should be stored in a subdirectory 
dataDir: test_data/varcall_test_data/raw
# outDir: output directory
outDir: test_data/varcall_test_data/output
# sampleFile: text file listing samples to be analysed
sampleFile: test_data/varcall_test_data/test_samples.txt
# forward and reverse fastq suffixes
fq_fwd: _R1.fq.gz
fq_rvr: _R2.fq.gz

# Align to reference
# Reference genome to call variants against
reference: test_data/varcall_test_data/LL6_1.fasta.gz

# The rest of the parameters don't have to be changed, copy and paste into your config file
# Preprocessing 
qc: yes # yes to perform preprocessing, no to skip

# BBMap paramerters
mink: 11
trimq: 14
mapq: 20
minlen: 45

merged: false # By default do not use merge for isolate genome assembly
fastqc: no # Run fastqc or not. Options: no, before, after, both. 


# Standard parameters. 
adapters: ../../../data/adapters/adapters.fa
phix: ../../../data/adapters/phix174_ill.ref.fa.gz

```

Example structure for the `dataDir`:

```
.
└── dataDir/
    ├── Sample1/
    │   ├── Sample1_ABCD_R1.fq.gz
    │   └── Sample1_EFGH_R2.fq.gz
    └── Sample2/
        ├── Sample2_ABCD_R1.fq.gz
        └── Sample2_EFGH_R2.fq.gz
```

In which case, `sampleFile` would look like this:

```
Sample1
Sample2
```

#### 3. Running variant calling pipeline. (Recommened to run with `--dry` option first)

```bash
nccrPipe isolate -c <full/path/to/your_config.yaml> -m call_variants  
```


#### 4. Running genome assembly pipeline.

```bash
nccrPipe isolate -c <full/path/to/your_config.yaml> -m assemble
```


#### Troubleshooting

Log files, stdout and stderr files for each step of the pipeline can be found in `outDir/logs/`


----------------------------------------------------------
# UNDER CONSTRUCTION 

### Analysis Options:

- preprocess
- merge_fastq
- assemble
- quastCheck
- plasmid
- annotate
- nucmer
- runANI
- phylogeny
- findAb
- pileup
- align
- align_with_ref
- call_vars
- call_vars2 -> use this one
- call_vars3
- markdup
- type
- serotype
- anVar
- anVar2
- pangenome

All of these need to be tested and documented

## Example

Default data structure: 


By default, data will be in `data/raw` and the output directory `data/processed`

```
cd project
nccrPipe -a create_config -c code/configs/project_config.yaml

```


## Variant Calling

### General Steps:
1. Align reads to reference genome with BWA. 
2. Remove duplicates with GATK MarkDuplicates
3. Run `bcftools mpileup` + `bcftools call`
4. `bcftools filter` 

```
-g5    filter SNPs within 5 base pairs of an indel or other other variant type
-G10   filter clusters of indels separated by 10 or fewer base pairs allowing only one to pass
-e  exclude
QUAL<10 calls with quality score < 10
DP4[2]<10 || DP4[3]<10  calls with < 10 reads (forward and reverse) covering the variant
(DP4[2] + DP4[3])/sum(DP4) < 0.9 calls with allele frequence < 90 %
MQ<50 calls with average mapping quality < 50 

```  
- Add filter based on coverage? Regions with really high coverage generally contain lots of artifacts. 

### Filters:


### Annotation:

- Still a little complicated. Adding new genome to the snpEff database not integrated into the main pipeline.
- Created conda environment within the pipeline: `.snakemake/conda/`
- If want to add new genome, run
```
nccrPipe -c <config_file> -a addGenomeSnpEff
```

- Config file has to include path to the gbk (and ideally fasta) files, as well as the name of the genome. The chromosome names between the files have to match. 
- If that doesn't fail, can run

```
nccrPipe -c <config_file> -a anVar
```

a
