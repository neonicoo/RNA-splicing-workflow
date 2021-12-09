# RNA-splicing-workflow
## Project 3: Implement a bioinformatic workflow to detect differential RNA splicing

[![CircleCI](https://circleci.com/gh/IARCbioinfo/template-nf.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/template-nf)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/template-nf/)
[![DOI](https://zenodo.org/badge/94193130.svg)](https://zenodo.org/badge/latestdoi/94193130)

![Workflow representation](template-nf.png)

## Description
### Background

Many cancers are at least partly initiated by the expression of aberrant transcript isoforms. RNA-seq is often analyzed using known annotated isoforms but tumors can also express novel isoforms that can only detected using dedicated tools.

### Data

- simple fastq file, available at https://github.com/IARCbioinfo/data_test

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
  - [trim-galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  - [salmon](https://combine-lab.github.io/salmon/)
  - [R](https://www.r-project.org/)
  - [SUPPA2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1417-1)

You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.

## Workflow

1. Paired-end reads trimming with trim-galore
2. RNA-seq quantification with salmon
3. Splicing analysis with SUPPA2

# Files and Warnings (WIP)

- Index building with salmon needs to be performed on the **transcriptome** not the genome (`--ref` parameter, [transcriptome that can be used](https://github.com/comprna/SUPPA_supplementary_data/blob/master/annotation/hg19_EnsenmblGenes_sequence_ensenmbl.fasta.gz))
- R [script](https://github.com/comprna/SUPPA/blob/master/scripts/format_Ensembl_ids.R) to format transcripts ids (`--formatscript` parameter)
- Event calculation (`--annot` parameter) requires an unziped GTF file ([one that can be used from Ensembl](https://github.com/comprna/SUPPA_supplementary_data/blob/master/annotation/Homo_sapiens.GRCh37.75.formatted.gtf.gz))

## Input
  | Type      | Description     |
  |-----------|---------------|
  | input1    | ...... |
  | input2    | ...... |

  Specify the test files location

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --param1    |            xx | ...... |
| --param2    |            xx | ...... |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --param3   |            xx | ...... |
| --param4    |            xx | ...... |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |
| --flag2    |      .... |


## Usage
  ```
  ...
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | output1    | ...... |
  | output2    | ...... |
