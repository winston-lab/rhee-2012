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
        expand("qual_ctrl/fastqc/{sample}/{sample}-fastqc_raw.html", sample=list(SAMPLES.keys())+["unknown"]),
        expand("alignment/{sample}.bam", sample=SAMPLES)

rule make_barcode_fasta:
    output:
        "fastq/barcodes.fa"
    run:
        with open(output[0], "w") as out:
            for k,v in SAMPLES.items():
                out.write('>'+k+'\n^'+v["barcode"]+'\n')

#minimum length 5 sanitizes output for bowtie 
rule trim_and_demultiplex:
    input:
        fq = config["fastq"],
        barcodes = "fastq/barcodes.fa"
    params:
        trim_qual = config["cutadapt"]["trim_qual"]
    output:
        temp(expand("fastq/{sample}.fastq", sample=list(SAMPLES.keys())+["unknown"]))
    log: "logs/trim_and_demultiplex.log"
    shell: """
        cutadapt -c -g file:{input.barcodes} -q {params.trim_qual} -m 5 -o fastq/{{name}}.fastq --format=sra-fastq --length-tag 'length=' {input.fq}
        """

rule fastqc:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        html = "qual_ctrl/fastqc/{sample}/{sample}-fastqc_raw.html",
        folder = "qual_ctrl/fastqc/{sample}/{sample}-fastqc_raw.zip",
    threads: config["threads"]
    log: "logs/fastqc_raw/fastqc-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.sample}) &> {log}
        (fastqc -o qual_ctrl/fastqc/{wildcards.sample} --noextract -t {threads} {input}) &>> {log}
        """

#align with bowtie1 on colorspace index
rule bowtie_build:
    input:
        fasta = config["genome"]["fasta"]
    output:
        expand("{idx_path}/{basename}.{num}.ebwt", idx_path=config["bowtie"]["index-path"], basename=config["genome"]["name"], num=[1,2,3,4]),
        expand("{idx_path}/{basename}.rev.{num}.ebwt", idx_path=config["bowtie"]["index-path"], basename=config["genome"]["name"], num=[1,2])
    params:
        idx_path = config["bowtie"]["index-path"],
        prefix = config["genome"]["name"]
    log: "logs/bowtie_build.log"
    shell: """
        (bowtie-build --color {input.fasta} {params.idx_path}/{params.prefix}) &> {log}
        """

#using bowtie 1.1.2 because 1.2.1 has an error with colorspace reads
#annoyingly, 1.1.2 doesn't handle gzip input so we have to zip afterwards
rule align:
    input:
        expand("{idx_path}/{basename}.{num}.ebwt", idx_path=config["bowtie"]["index-path"], basename=config["genome"]["name"], num=[1,2,3,4]),
        expand("{idx_path}/{basename}.rev.{num}.ebwt", idx_path=config["bowtie"]["index-path"], basename=config["genome"]["name"], num=[1,2]),
        fastq = "fastq/{sample}.fastq"
    params:
        idx_path = config["bowtie"]["index-path"] + "/" + config["genome"]["name"],
        max_mismatch = config["bowtie"]["max_mismatch"]
    output:
        bam = "alignment/{sample}.bam",
        unaligned = temp("alignment/unaligned-{sample}.fastq")
    conda: "envs/bowtie.yaml"
    threads: config["threads"]
    log:
        "logs/align/align-{sample}.log"
    shell: """
        (bowtie --color -v {params.max_mismatch} --nomaqround --best -S -p {threads} --un {output.unaligned} {params.idx_path} {input.fastq} | samtools view -buh -q 50 - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {log}
        """

rule gzip_loose_fastq:
    input:
        "alignment/{sample}.bam",
        cleaned =  "fastq/{sample}.fastq",
        unaligned = "alignment/unaligned-{sample}.fastq"
    output:
        "fastq/{sample}.fastq.gz",
        "alignment/unaligned-{sample}.fastq.gz"
    shell: """
        pigz -fk {input.cleaned}
        pigz -fk {input.unaligned}
        """











