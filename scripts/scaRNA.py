from genericpath import exists
import json
import os
import re
import sys
from pathlib import Path
from subprocess import call, check_output
from subprocess import run, PIPE
import logging
import argparse

import mappy as mp
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("scaRNA")

project_root = Path(__file__).resolve().parents[1]


data_dirs = {}
for folder in ["trim", "alignment", "rRNA", "bed"]:
    data_dirs[folder] = project_root / "data" / "interim" / folder
    if not data_dirs[folder].exists():
        logger.info(f"Creating {data_dirs[folder]}")
        data_dirs[folder].mkdir(parents=True)

for folder in ["bigwig", "bedgraph", "json", "peak"]:
    data_dirs[folder] = project_root / "data" / "output" / folder
    if not data_dirs[folder].exists():
        logger.info(f"Creating {data_dirs[folder]}")
        data_dirs[folder].mkdir(parents=True)


def run_cmd(cmd):
    logger.info(f"Running command: {cmd}")
    res = run(cmd, shell=True, check=True, stdout=PIPE, stderr=PIPE)
    if res.returncode != 0:
        logger.error(f"Error running command: {cmd}")
        logger.error(f"Error message: {res.stderr.decode()}")
        sys.exit(1)


def trim(sample, fq1, fq2, three_primer="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", primer_match_len=10):
    logger.info(f"Trimming small RNA Kit primer for {sample}")
    three_primer = three_primer[:primer_match_len]
    total = 0
    with_primer, without_primer = 0, 0
    yes_insert, no_insert = 0, 0
    trim_fq1_name = data_dirs["trim"] / f"{sample}_trim_1.fq"
    trim_fq2_name = data_dirs["trim"] / f"{sample}_trim_2.fq"
    trim_fq1 = open(trim_fq1_name, "w")
    trim_fq2 = open(trim_fq2_name, "w")
    for read1, read2 in zip(
        mp.fastx_read(fq1, read_comment=True),
        mp.fastx_read(fq2, read_comment=True),
    ):
        total += 1
        read_name, comment = read1[0], read1[-1][1:]
        seq1, qual1 = read1[1], read1[2]
        seq2, qual2 = read2[1], read2[2]
        if three_primer in seq1:
            with_primer += 1
            if seq1.startswith(three_primer):
                no_insert += 1
                continue
            else:
                yes_insert += 1
                insert_len = len(seq1.split(three_primer)[0])
        else:
            without_primer += 1
            insert_len = 150
        trim_fq1.write(f"@{read_name} 1{comment}\n{seq1[:insert_len]}\n+\n{qual1[:insert_len]}\n")
        trim_fq2.write(f"@{read_name} 2{comment}\n{seq2[:insert_len]}\n+\n{qual2[:insert_len]}\n")
        if total % 1000000 == 0:
            logger.info(
                f"Processed {total} reads, with primer: {with_primer}, without primer: {without_primer}"
            )
    trim_fq1.close()
    trim_fq2.close()

    trim_count = {
        "total": total,
        "with_primer": with_primer,
        "without_primer": without_primer,
        "yes_insert": yes_insert,
        "no_insert": no_insert,
    }
    logger.info(f"Trimming complete for {sample}")
    logger.info(f"Trim results: {trim_count}")
    with open(data_dirs["json"] / f"{sample}.trim.json", "w") as f:
        json.dump(trim_count, f, indent=4)
    logger.info(f"Compressing trimmed files for {sample}")
    run_cmd(f"gzip -f {trim_fq1_name}")
    run_cmd(f"gzip -f {trim_fq2_name}")


def filter_rRNA(sample):
    rRNA_dir = data_dirs["rRNA"]
    if not os.path.exists(f"{ref_dir}/rRNA_databases_v4/smr_v4.3_default_db.fasta"):
        os.makedirs(f"{ref_dir}/rRNA_databases_v4", exist_ok=True)
        logger.info("rRNA databases not found")
        url = "https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz"
        logger.info(f"Downloading rRNA databases from {url}")
        run_cmd(f"wget {url} -O {ref_dir}/database.tar.gz")
        run_cmd(f"tar -xvf {ref_dir}/database.tar.gz -C {ref_dir}/rRNA_databases_v4")
        if not os.path.exists(f"{ref_dir}/rRNA_databases_v4/smr_v4.3_default_db.fasta"):
            logger.error("failed to download rRNA databases")
            sys.exit(1)
    trim_fq1_name = data_dirs["trim"] / f"{sample}_trim_1.fq.gz"
    trim_fq2_name = data_dirs["trim"] / f"{sample}_trim_2.fq.gz"
    kvdb = rRNA_dir / "kvdb" / sample
    if kvdb.exists():
        logger.info(f"Removing existing kvdb {kvdb}")
        run_cmd(f"rm -rf {kvdb}")
    run_cmd(f"mkdir -p {kvdb}")
    cmd = f"""sortmerna \
    --threads 48 \
    --kvdb {kvdb} \
    --ref {ref_dir}/rRNA_databases_v4/smr_v4.3_default_db.fasta \
    --workdir {ref_dir}/rRNA_databases_v4 \
    --reads {trim_fq1_name} \
    --reads {trim_fq2_name} \
    --aligned "{sample}/N_rRNA" --other "{sample}/N_clean" --paired_in --fastx --out2"""
    run_cmd(cmd)
    call(f"cp ./{sample} {rRNA_dir}", shell=True)
    run_cmd(f"rm -rf {sample}")


def align(sample, fq1, fq2):
    bowtie2_index = f"{ref_dir}/Bowtie2Index/genome"
    if not os.path.exists(f"{bowtie2_index}.1.bt2"):
        if genome == 'hg19':
            species = 'Homo_sapiens'
        elif genome == 'mm10':
            species = 'Mus_musculus'
        else:
            logger.error('Unsupported genome')
            sys.exit(1)
        logger.info(f"{bowtie2_index} not found")
        logger.info("Downloading from iGenomes")
        cmd = f"""aws s3 \
        --no-sign-request \
        --region eu-west-1 sync \
        s3://ngi-igenomes/igenomes/{species}/UCSC/{genome}/Sequence/Bowtie2Index/ \
        {ref_dir}/Bowtie2Index/
        """
        run_cmd(cmd)
        if not os.path.exists(f"{bowtie2_index}.1.bt2"):
            logger.error("Failed to download bowtie2 index")
            sys.exit(1)
    logger.info("Align with bowtie2")
    cmd = f"""bowtie2 -p 10 --very-sensitive \
        -x {bowtie2_index} -1 {fq1} -2 {fq2} \
        |samtools view -bS - > \
        {data_dirs["alignment"]}/{sample}.bam"""
    run_cmd(cmd)
    logger.info("Sort bam file by read names.")
    cmd = f"""samtools sort -n \
    -@ 8 -m 4G \
    -O bam -o {data_dirs["alignment"]}/{sample}.sorted.bam \
    {data_dirs["alignment"]}/{sample}.bam"""
    run_cmd(cmd)
    logger.info("Sort bam file by coordinates.")
    cmd = f"""samtools sort \
        -@ 8 -m 4G -O bam \
        -o {data_dirs["alignment"]}/{sample}.hg19.bam \
        {data_dirs["alignment"]}/{sample}.bam"""
    run_cmd(cmd)
    logger.info("index sorted bam file")
    cmd = f"""samtools index {data_dirs["alignment"]}/{sample}.hg19.bam"""
    run_cmd(cmd)


def prepare_beds():
    gtf_file = f'{ref_dir}/gencode.gtf.gz'
    if not os.path.exists(gtf_file):
        logger.info("Downloading gencode gtf file")
        if genome == 'hg19':
            url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
        elif genome == 'mm10':
            url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz"
        else:
            logger.error('Unsupported genome build')
            sys.exit(1)
        cmd = f"""wget {url} -O {gtf_file}"""
        run_cmd(cmd)
    biotypes = ["miRNA", "snoRNA", "snRNA", "mRNA", "mRNA_tss"]
    logger.info('Generating bed files for {}'.format(biotypes))
    for biotype in biotypes[:-1]:
        call(
            f'zgrep "{biotype}" {gtf_file} > {ref_dir}/{biotype}.bed',
            shell=True,
        )
        a = pd.read_csv(f"{ref_dir}/{biotype}.bed", sep="\t", header=None)
        a = a[a[2] == "transcript"]
        a[9] = 0
        a[[0, 3, 4, 5, 9, 6]].drop_duplicates().to_csv(
            f"{ref_dir}/{biotype}.bed", sep="\t", index=False, header=False
        )
    all_mRNA = pd.read_csv(gtf_file, sep="\t", header=None, comment="#")
    all_mRNA = all_mRNA[all_mRNA[2] == "transcript"]
    all_mRNA[9] = 0
    all_mRNA = all_mRNA[[0, 3, 4, 5, 9, 6]].drop_duplicates()
    all_mRNA.to_csv(f"{ref_dir}/mRNA.bed", sep="\t", index=False, header=False)
    all_mRNA = pd.read_csv(f"{ref_dir}/mRNA.bed", sep="\t", header=None)
    all_mRNA["tss"] = all_mRNA[1].where(all_mRNA[5] == "+", all_mRNA[2])
    all_mRNA[1] = all_mRNA["tss"] - 100
    all_mRNA[2] = all_mRNA["tss"] + 100
    all_mRNA[all_mRNA[1] > 0][[0, 1, 2, 3, 4, 5]].to_csv(
        f"{ref_dir}/mRNA_tss.bed", sep="\t", index=False, header=False
    )


def tobed(sample):
    logger.info("Convert bam to bedpe")
    cmd = f"""bedtools bamtobed \
            -bedpe -mate1 \
            -i {data_dirs["alignment"]}/{sample}.sorted.bam | \
            grep chr | \
            awk '{{if ($1==$4) print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$9}}' \
            > {data_dirs["alignment"]}/{sample}_bedpe_fragments.txt"""
    run_cmd(cmd)

    df = pd.read_csv(
        f'{data_dirs["alignment"]}/{sample}_bedpe_fragments.txt', sep="\t", header=None
    )
    plus = df[df[5] == "+"][[0, 1, 4]]
    plus[10] = 0
    plus[11] = 0
    plus[12] = "+"
    plus.columns = range(6)

    minus = df[df[5] == "-"][[0, 3, 2]]
    minus[10] = 0
    minus[11] = 0
    minus[12] = "-"
    minus.columns = range(6)
    df = pd.concat([minus, plus], ignore_index=True)
    df = df[(df[1] < df[2]) & ((df[2] - df[1]) <= 1000)]
    df = df.sort_values([0, 1])

    df.to_csv(
        f'{data_dirs["bed"]}/{sample}_frags.sort.bed',
        sep="\t",
        index=False,
        header=False,
    )


def intersect(sample):
    biotypes = ["miRNA", "snoRNA", "snRNA", "mRNA", "mRNA_tss"]
    for i, biotype in enumerate(biotypes[:-2]):
        if i == 0:
            abed = f'{data_dirs["bed"]}/{sample}_frags.sort.bed'
        else:
            abed = f'{data_dirs["bed"]}/{sample}_frags.sort.rm{biotypes[i-1]}.bed'
        logger.info(f"Remove {biotype}")
        cmd = f"""bedtools intersect \
            -a {abed} \
            -b {ref_dir}/{biotype}.bed \
            -u -names -s > \
            {data_dirs["bed"]}/{sample}_frags.sort.{biotype}.bed"""
        run_cmd(cmd)
        cmd = f"""bedtools intersect \
            -a {abed} \
            -b {ref_dir}/{biotype}.bed \
            -v -s > \
            {data_dirs["bed"]}/{sample}_frags.sort.rm{biotype}.bed"""
        run_cmd(cmd)

    for biotype in biotypes[-2:]:
        logger.info(f'Extract {biotype}')
        cmd = f"""bedtools intersect \
            -a {data_dirs["bed"]}/{sample}_frags.sort.rm{biotypes[-3]}.bed \
            -b {ref_dir}/{biotype}.bed \
            -u -s > \
            {data_dirs["bed"]}/{sample}_frags.sort.{biotype}.bed"""
        run_cmd(cmd)


def tss_tes(sample):
    bg_dir = data_dirs["bedgraph"]
    bw_dir = data_dirs["bigwig"]
    genome_file = f'{ref_dir}/genome.file'
    if not os.path.exists(genome_file):
        logger.info("Generate genome file from fa.fai")
        cmd = f"""cut -f1,2 {ref_dir}/Bowtie2Index/genome.fa.fai > {genome_file}"""
        run_cmd(cmd)
    bedtools_i = f'bedtools genomecov -i {data_dirs["bed"]}/{sample}_frags.sort.mRNA.bed -g {genome_file}'
    logger.info('Make strand specific bedgraph files for TSS and TES regions.')
    cmd = f'{bedtools_i} -strand + -5 -bg > {bg_dir}/{sample}.fwd.tss.bedgraph',
    run_cmd(cmd)
    cmd = f'{bedtools_i} -strand + -3 -bg > {bg_dir}/{sample}.fwd.tes.bedgraph',
    run_cmd(cmd)
    cmd = f'{bedtools_i} -strand - -5 -bg > {bg_dir}/{sample}.rev.tss.bedgraph',
    run_cmd(cmd)
    cmd = f'{bedtools_i} -strand - -3 -bg > {bg_dir}/{sample}.rev.tes.bedgraph',
    run_cmd(cmd)

    bgtobw = "bedGraphToBigWig"
    logger.info('Convert bedgraph to bigwig')
    cmd = f'{bgtobw} {bg_dir}/{sample}.fwd.tss.bedgraph {genome_file} {bw_dir}/{sample}.fwd.tss.bw'
    run_cmd(cmd)
    cmd = f'{bgtobw} {bg_dir}/{sample}.fwd.tes.bedgraph {genome_file} {bw_dir}/{sample}.fwd.tes.bw'
    run_cmd(cmd)
    cmd = f'{bgtobw} {bg_dir}/{sample}.rev.tss.bedgraph {genome_file} {bw_dir}/{sample}.rev.tss.bw'
    run_cmd(cmd)
    cmd = f'{bgtobw} {bg_dir}/{sample}.rev.tes.bedgraph {genome_file} {bw_dir}/{sample}.rev.tes.bw'
    run_cmd(cmd)


def tsspeak(sample):
    logger.info("Call tss peak")
    cmd = f"""makeTagDirectory \
        {data_dirs["peak"]}/{sample} \
        {data_dirs["alignment"]}/{sample}.hg19.bam \
        -genome {ref_dir}/Bowtie2Index/genome.fa \
        -checkGC -fragLength 100"""
    run_cmd(cmd)
    cmd = f"""findcsRNATSS.pl \
        {data_dirs["peak"]}/{sample} \
        -o {data_dirs["peak"]}/{sample}.hg19.peak \
        -gtf {ref_dir}/gencode.gtf.gz \
        -genome {ref_dir}/Bowtie2Index/genome.fa"""
    run_cmd(cmd)


def summary(sample):
    flagstat = check_output(
        f'samtools flagstat {data_dirs["alignment"]}/{sample}.sorted.bam', shell=True
    )
    total_read = int(re.findall(r"(\d+) \+ 0 read1", flagstat.decode())[0])
    with open(f'{data_dirs["rRNA"]}/{sample}/N_rRNA.log') as f:
        rrna = f.readlines()
    rrna = int(re.findall(r"Total reads passing E-value threshold = (\d+) ", "".join(rrna))[0]) // 2

    def rows_count(file):
        rows = check_output(f"wc -l {file}", shell=True).decode()
        rows = int(rows.split()[0])
        return rows

    total_pairs = rows_count(f'{data_dirs["bed"]}/{sample}_frags.sort.bed')
    rmmirna = rows_count(f'{data_dirs["bed"]}/{sample}_frags.sort.rmmiRNA.bed')
    rmsnorna = rows_count(f'{data_dirs["bed"]}/{sample}_frags.sort.rmsnoRNA.bed')
    rmsnrna = rows_count(f'{data_dirs["bed"]}/{sample}_frags.sort.rmsnRNA.bed')
    mirna = total_pairs - rmmirna
    snorna = rmmirna - rmsnorna
    snrna = rmsnorna - rmsnrna
    mrna = rows_count(f'{data_dirs["bed"]}/{sample}_frags.sort.mRNA.bed')
    tss = rows_count(f'{data_dirs["bed"]}/{sample}_frags.sort.mRNA_tss.bed')
    mrna = mrna - tss
    other = total_pairs - sum([mirna, snorna, snrna, mrna, tss])
    map_rrna = {
        "unmapped": total_read - total_pairs,
        "rRNA": rrna,
        "miRNA": mirna,
        "snoRNA": snorna,
        "snRNA": snrna,
        "other_mapped": other,
        "mRNA": mrna,
        "mRNA_tss": tss,
    }
    with open(f'{data_dirs["json"]}/{sample}.summary.json', "w") as f:
        json.dump(map_rrna, f, indent=4)


def parse_args():
    parser = argparse.ArgumentParser(description="scaRNA-seq analysis")
    parser.add_argument("--sample", type=str, help="Sample name")
    parser.add_argument("--fq1", type=str, help="Read1 fastq file")
    parser.add_argument("--fq2", type=str, help="Read2 fastq file")
    parser.add_argument(
        "--genome", type=str, default="hg19", help="Genome version", choices=["hg19", "mm10"]
    )
    parser.add_argument("--three_primer", type=str, default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC")
    return parser.parse_args()


def main():
    args = parse_args()
    sample = args.sample
    fq1 = args.fq1
    fq2 = args.fq2
    global genome
    genome = args.genome
    global ref_dir
    ref_dir = project_root / "data" / "ref" / genome
    ref_dir.mkdir(parents=True, exist_ok=True)
    trim(sample, fq1, fq2)
    filter_rRNA(sample)
    align(
        sample,
        f'{data_dirs["rRNA"]}/{sample}/N_clean_fwd.fq.gz',
        f'{data_dirs["rRNA"]}/{sample}/N_clean_rev.fq.gz',
    )
    tobed(sample)
    if not os.path.exists(f"{ref_dir}/mRNA_tss.bed"):
        prepare_beds()
    intersect(sample)
    tss_tes(sample)
    summary(sample)
    tsspeak(sample)


if __name__ == "__main__":
    main()
