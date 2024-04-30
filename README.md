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
git clone https://github.com/your-repository/scaRNA-seq-tool.git
cd scaRNA-seq-tool
```
create conda environment:
```
conda env create -f environment.yml
```
## Usage

To run the tool, you need to provide sample names, paths to the read files (both forward and reverse), and specify the genome version. The command uses the following structure:

```
python script/scaRNA.py \
--sample [sample_name] \
--fq1 [path_to_read1.fq] \
--fq2 [path_to_read2.fq] \
--genome [hg19/mm10] \
--ref_dir [path_to_reference_directory]
```

Optional Arguments
- --three_primer: A string representing the adapter sequence for three-prime trimming (default is specific to the kit used).
Example Command
python script.py --sample sample1 --fq1 path/to/sample1_R1.fq --fq2 path/to/sample1_R2.fq --genome hg19 --ref_dir path/to/ref

 
Output

The tool generates multiple output files organized in specific directories under the data folder. Key outputs include:

Trimmed reads: Compressed files containing reads after trimming.
Filtered reads: Reads after rRNA filtering.
Aligned reads: BAM files after genomic alignment.
Annotation files: BED files corresponding to different genomic features.
Summary files: JSON files containing summary statistics of processing steps.
Contributing

Contributions to this project are welcome! Please fork the repository and submit pull requests with any enhancements or bug fixes.

License

This project is licensed under the MIT License - see the LICENSE file for details.

Contact Information

For help or feedback, please file an issue in the GitHub repository.

This README provides you with detailed instructions on how to use the scaRNA-seq analysis tool and outlines the features and workflow of the script comprehensively.