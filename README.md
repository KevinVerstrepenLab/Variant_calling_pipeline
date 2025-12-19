# Variant Calling Pipeline (Snakemake)

This repository contains a Snakemake workflow for **FASTQ → BAM → GVCF → joint genotyping → filtered VCFs** using **GATK4**, **BWA-MEM**, **fastp**, and optional downstream filtering.  
It is designed for reproducible variant calling of multiple samples on an HPC cluster or workstation.

## Pipeline Overview
1. Quality trimming with fastp
2. Alignment to reference genome using bwa-mem2
3. Sorting & duplicate marking
4. GVCF creation per sample (GATK HaplotypeCaller)
5. GenomicsDBImport (joint calling)
6. Joint genotyping
7. Variant filtering
8. Final VCFs with missingness filtering

## Repository Structure
- config/ — config.yaml and samples.tsv  
- workflow/ — Snakefile and rule modules  
- workflow/rules/ — individual .smk rule files  
- cmd.unix — optional cluster command example

## Sample Sheet Format
sample<TAB>R1<TAB>R2  
S1   S1_R1.fq.gz   S1_R2.fq.gz

## Running the Pipeline
```
snakemake --use-conda --cores 56 --printshellcmds --latency-wait 30           --keep-going --rerun-incomplete           --configfile config/config.yaml
```

## Outputs
- BAM files
- GVCFs
- Joint-called VCF
- Filtered SNP and INDEL VCFs

## Troubleshooting
- `snakemake -n` to dry-run  
- `snakemake --debug-dag` to visualize  
- `snakemake --rulegraph | dot -Tpng > dag.png` for workflow graph


## Maintainer
Verstrepen Lab — KU Leuven
michael.abrouk@kuleuven.be
