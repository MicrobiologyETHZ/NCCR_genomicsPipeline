# NCCR Genomics Pipeline

**Under Development**


Pipeline for analysis of isolate genomes. 

### Steps Currently Included:

- QC
- Isolate Genome Assembly 
- Gene calling and annotation
- Mapping back to assembly or to reference genome
- Variant calling and annotation
- ANI calculation
- Phylogenetics with PhyloPhlAn
- Pangenome analysis with panX


### Installation

The software depends on `conda`, `Python >= 3.8` and `snakemake == 5.22`. 

```
git clone https://github.com/MicrobiologyETHZ/NCCR_genomicsPipeline
cd code/package
conda env create -f nccrPipe_environment.yml
conda activate nccrPipe
pip install -e . 

```

### Running

```
nccrPipe -c <config_file> -a <analysis> [-np]

    -c Path to config file (default=test config file)
    -a Analysis to run, ex. assemble, call_vars (under development)
    -np Dry run. Will run snakemake with -np flag
```

To start a new project, run 

```
nccrPipe -c <config_file> -a create_config
```
This will copy the default version of config to `<config_file>`. Right now this file needs to be manually edited to input relevant info about the project.
After the `config_file` is setup, run the analysis. For analysis options, see below. 

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