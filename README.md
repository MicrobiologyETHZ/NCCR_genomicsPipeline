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

## Example

Default data structure: 

```
|-- project/
    |--code/
    |  |--configs/
    |     |--config.yaml
    |     |--samples.txt 
    | 
    |--data/
       |--raw/
       |--processed/     

```
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