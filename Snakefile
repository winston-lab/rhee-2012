#!/usr/bin/env python

configfile: "config.yaml"

SAMPLES=config["samples"]

localrules: all,
    make_barcode_fasta

rule all:
    input:
        #FastQC
        "fastq/barcodes.fa",
        expand("fastq/{sample}.fastq.gz", sample=list(SAMPLES.keys())+["unknown"]),
        expand("qual_ctrl/fastqc/cleaned/{sample}/{sample}-fastqc_raw.html", sample=list(SAMPLES.keys())+["unknown"])


rule make_barcode_fasta:
    output:
        "fastq/barcodes.fa"
    run:
        with open(output[0], "w") as out:
            for k,v in SAMPLES.items():
                out.write('>'+k+'\n^'+v["barcode"]+'\n')

rule trim_and_demultiplex:
    input:
        fq = config["fastq"],
        barcodes = "fastq/barcodes.fa"
    params:
        trim_qual = config["cutadapt"]["trim_qual"]
    output:
        expand("fastq/{sample}.fastq.gz", sample=list(SAMPLES.keys())+["unknown"])
    log: "logs/trim_and_demultiplex.log"
    shell: """
        cutadapt -c -g file:{input.barcodes} -q {params.trim_qual} -o fastq/{{name}}.fastq.gz --format=sra-fastq --length-tag 'length=' {input.fq}
        """

rule fastqc:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        html = "qual_ctrl/fastqc/cleaned/{sample}/{sample}-fastqc_raw.html",
        folder = "qual_ctrl/fastqc/cleaned/{sample}/{sample}-fastqc_raw.zip",
    threads: config["threads"]
    log: "logs/fastqc_raw/fastqc_raw-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/cleaned/{wildcards.sample}) &> {log}
        (fastqc -o qual_ctrl/fastqc/cleaned/{wildcards.sample} --noextract -t {threads} {input}) &>> {log}
        """
