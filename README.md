# scaRNA-seq Analysis Pipeline
## Overview

This tool is designed to process small RNA sequencing data in a systematic fashion. The workflow includes trimming of adaptors, filtering ribosomal RNA (rRNA), alignment to a reference genome, and post-processing steps such as generation of bed files and peak calling. It outputs various intermediate files, such as trimmed reads, aligned bam files, and summary statistics in JSON format.

## Features
- Trimming of Small RNA Kit Primer: Customizable trimming of adaptors based on the provided primer sequence.
- rRNA Filtering: Uses SortMeRNA for filtering rRNA from small RNA-seq data.
- Alignment: Utilizes Bowtie2 for alignment to a reference genome.
- Converter to BED: Converts BAM files to BED format for easier manipulation of genomic intervals.
- Remove lncRNA: remove miRNA, snoRNA, and snRNA from the annotation file
- TSS and TES analysis: Generates specific files for Transcription Start Sites (TSS) and Transcription End Sites (TES).

## Installation

Clone the repository to your local machine:
```
git clone https://github.com/Jianhua-Wang/scaRNAseq.git
cd scaRNAseq
```
create conda environment:
```
conda env create -f environment.yml
conda activate scarna
```
## Usage

To run the tool, you need to provide sample names, paths to the read files (both forward and reverse), and specify the genome version. The command uses the following structure:

```
python scripts/scaRNA.py \
--sample [sample_name] \
--fq1 [path_to_read1.fq] \
--fq2 [path_to_read2.fq] \
--genome [hg19/mm10] \
--three_primer [adapter_sequence]
```
The idx of `sortmerna` and `bowtie2` will be downloaded automatically. The idx files are stored in `data/ref` directory.

scaRNA-seq involves small RNA library construction, which is a directional library. If the kit works fine, the cDNA fragment will be consist of 5' primer, insert cDNA, and 3' primer. So, the reverse read will start with the 3' primer when there is no insert. That's way we need to specify the 3' primer sequence to trim the 3' primer sequence from the reverse read.

## Example
```
python scripts/scaRNA.py \
--sample test \
--fq1 data/raw/test_1.fq.gz \
--fq2 data/raw/test_2.fq.gz \
--genome hg19 \
--three_primer AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```
In this case, we used a example dataset in the `data/raw` directory. The forward read is `test_1.fq.gz` and the reverse read is `test_2.fq.gz`. The genome version is `hg19` and the 3' primer sequence is `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`.

You can generate the QC plots using `./scripts/scaRNA_summary.ipynb`.

## Output

The tool generates multiple output files organized in specific directories under the data folder. Key outputs include:

- Trimmed reads: Compressed files containing reads after trimming.
- Filtered reads: Reads after rRNA filtering.
- Aligned reads: BAM files after genomic alignment.
- Annotation files: BED files corresponding to different genomic features.
- Summary files: JSON files containing summary statistics of processing steps.