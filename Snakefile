#!/usr/bin/env python

configfile: "config.yaml"

SAMPLES=config["samples"]
NEXUS_SAMPLES=config["nexus_samples"]

localrules: all,
    make_barcode_fasta, make_stranded_genome

rule all:
    input:
        #FastQC
        "csfasta/barcodes.fa",
        expand("qual_ctrl/fastqc/{sample}/{sample}_fastqc.html", sample=list(SAMPLES.keys())+["unknown"]),
        expand("csfasta/{sample}.csfasta.gz", sample=list(SAMPLES.keys())+["unknown"]),
        expand("alignment/{sample}.bam", sample=list(SAMPLES.keys())+["unknown"]),
        expand("coverage/{norm}/{sample}-chipexo-{norm}-{strand}.{fmt}", norm=["counts","libsizenorm"], strand=["plus", "minus", "SENSE", "ANTISENSE"], sample=list(SAMPLES.keys())+["unknown"], fmt=["bedgraph", "bw"]),
        # expand("datavis/{annotation}/{annotation}-{sample}-libsizenorm-{strand}.tsv", annotation=config["annotations"], sample=list(SAMPLES.keys())+list(NEXUS_SAMPLES.keys()), strand=["SENSE", "ANTISENSE"]),
        expand("datavis/{annotation}/allsamples-{annotation}-libsizenorm-{strand}.tsv.gz", annotation=config["annotations"], strand=["SENSE", "ANTISENSE"]),
        expand("datavis/{annotation}/{annotation}-libsizenorm-exo-v-nexus-metagene.svg", annotation=config["annotations"])

rule make_barcode_fasta:
    output:
        "csfasta/barcodes.fa"
    run:
        with open(output[0], "w") as out:
            for k,v in SAMPLES.items():
                out.write('>'+k+'\n^'+v["barcode"]+'\n')

#minimum length 5 sanitizes output for bowtie 
#output is uncompressed because I don't feel like learning how to open gzip files in perl right now...
rule trim_and_demultiplex:
    input:
        csfasta = config["csfasta"],
        qual = config["qual"],
        barcodes = "csfasta/barcodes.fa"
    params:
        trim_qual = config["cutadapt"]["trim_qual"]
    output:
        temp(expand("csfasta/{sample}.csfastq", sample=list(SAMPLES.keys())+["unknown"]))
    log: "logs/trim_and_demultiplex.log"
    shell: """
        cutadapt -c -g file:{input.barcodes} -q {params.trim_qual} -m 5 -o csfasta/{{name}}.csfastq --length-tag 'length=' {input.csfasta} {input.qual}
        """

rule fastqc:
    input:
        "csfasta/{sample}.csfastq"
    output:
        html = "qual_ctrl/fastqc/{sample}/{sample}_fastqc.html",
        folder = "qual_ctrl/fastqc/{sample}/{sample}_fastqc.zip",
    threads: config["threads"]
    log: "logs/fastqc/fastqc-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.sample}) &> {log}
        (fastqc -o qual_ctrl/fastqc/{wildcards.sample} --noextract -t {threads} {input}) &>> {log}
        """

#cutadapt outputs in fastq format, but bowtie throws an error when I try to align using the fastq (either "Reads file contained a pattern with more than 1024 quality values" using bowtie 1.2.1.1, or "Too few quality values for read: SRR..." with bowtie 1.1.2) 
# so, convert back to csfasta/qual using a script I found on seqanswers by Nils Homer
rule csfastq_to_csfasta:
    input:
        "csfasta/{sample}.csfastq"
    output:
        temp("csfasta/{sample}.csfasta"),
        temp("csfasta/{sample}_QV.qual")
    shell: """
        perl scripts/csfastq_to_cs.pl {input}
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

rule align:
    input:
        expand("{idx_path}/{basename}.{num}.ebwt", idx_path=config["bowtie"]["index-path"], basename=config["genome"]["name"], num=[1,2,3,4]),
        expand("{idx_path}/{basename}.rev.{num}.ebwt", idx_path=config["bowtie"]["index-path"], basename=config["genome"]["name"], num=[1,2]),
        csfasta = "csfasta/{sample}.csfasta",
        qual = "csfasta/{sample}_QV.qual"
    params:
        idx_path = config["bowtie"]["index-path"] + "/" + config["genome"]["name"],
        max_mismatch = config["bowtie"]["max_mismatch"]
    output:
        bam = "alignment/{sample}.bam",
        unaligned = temp("alignment/unaligned-{sample}.fastq")
    # conda: "envs/bowtie.yaml"
    threads: config["threads"]
    log:
        "logs/align/align-{sample}.log"
    shell: """
        (bowtie --color -v {params.max_mismatch} --nomaqround --best -S -p {threads} --un {output.unaligned} {params.idx_path} -f {input.csfasta} --quals {input.qual} | samtools view -buh -q 50 - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {log}
        """

rule gzip_loose_files:
    input:
        "alignment/{sample}.bam",
        "qual_ctrl/fastqc/{sample}/{sample}_fastqc.html",
        "qual_ctrl/fastqc/{sample}/{sample}_fastqc.zip",
        csfasta = "csfasta/{sample}.csfasta",
        qual = "csfasta/{sample}_QV.qual",
        unaligned = "alignment/unaligned-{sample}.fastq"
    output:
        "csfasta/{sample}.csfasta.gz",
        "csfasta/{sample}_QV.qual.gz",
        "alignment/unaligned-{sample}.fastq.gz"
    shell: """
        pigz -f {input.csfasta}
        pigz -f {input.qual}
        pigz -f {input.unaligned}
        """

rule get_coverage:
    input:
        "alignment/{sample}.bam",
    output:
        plmin = "coverage/counts/{sample}-chipexo-counts-combined.bedgraph",
        plus = "coverage/counts/{sample}-chipexo-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-chipexo-counts-minus.bedgraph",
    log: "logs/get_coverage/get_coverage-{sample}.log"
    shell: """
        (genomeCoverageBed -bga -5 -ibam {input} | sort -k1,1 -k2,2n > {output.plmin}) &> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | sort -k1,1 -k2,2n > {output.minus}) &>> {log}
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | sort -k1,1 -k2,2n > {output.plus}) &>> {log}
        """

rule normalize:
    input:
        plmin = "coverage/counts/{sample}-chipexo-counts-combined.bedgraph",
        plus = "coverage/counts/{sample}-chipexo-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-chipexo-counts-minus.bedgraph",
    output:
        plus = "coverage/libsizenorm/{sample}-chipexo-libsizenorm-plus.bedgraph",
        minus = "coverage/libsizenorm/{sample}-chipexo-libsizenorm-minus.bedgraph",
    log: "logs/normalize/normalize-{sample}.log"
    shell: """
        (bash scripts/libsizenorm.sh {input.plmin} {input.plus} 1 > {output.plus}) &> {log}
        (bash scripts/libsizenorm.sh {input.plmin} {input.minus} 1 > {output.minus}) &>> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = "coverage/{norm}/{sample}-chipexo-{norm}-plus.bedgraph",
        minus = "coverage/{norm}/{sample}-chipexo-{norm}-minus.bedgraph",
    output:
        sense = "coverage/{norm}/{sample}-chipexo-{norm}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}-chipexo-{norm}-ANTISENSE.bedgraph"
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus} > {output.antisense}) &>> {log}
        """

rule make_stranded_genome:
    input:
        exp = config["genome"]["chrsizes"],
    output:
        exp = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
    log: "logs/make_stranded_genome.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.exp} > {output.exp}) &> {log}
        """

rule bg_to_bw:
    input:
        bedgraph = "coverage/{norm}/{sample}-chipexo-{norm}-{strand}.bedgraph",
        chrsizes = lambda wildcards: config["genome"]["chrsizes"] if wildcards.strand in ["plus","minus"] else os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
    output:
        "coverage/{norm}/{sample}-chipexo-{norm}-{strand}.bw",
    log: "logs/bg_to_bw/bg_to_bw-{sample}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
        """

rule make_stranded_annotations:
    input:
        lambda wildcards : config["annotations"][wildcards.annotation]["path"]
    output:
        "{annopath}/stranded/{annotation}-STRANDED.{ext}"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

rule deeptools_matrix:
    input:
        annotation = lambda wildcards: config["annotations"][wildcards.annotation]["path"] if wildcards.strand=="qfrags" else os.path.dirname(config["annotations"][wildcards.annotation]["path"]) + "/stranded/" + wildcards.annotation + "-STRANDED" + os.path.splitext(config["annotations"][wildcards.annotation]["path"])[1],
        bw = lambda wildcards: "coverage/libsizenorm/" + wildcards.sample +"-chipexo-libsizenorm-"+ wildcards.strand + ".bw" if wildcards.sample in SAMPLES else config["nexus_coverage_path"] + wildcards.sample +"-tfiib-chipnexus-libsizenorm-" + wildcards.strand + ".bw"
    output:
        dtfile = temp("datavis/{annotation}/{annotation}-{sample}-libsizenorm-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/{annotation}-{sample}-libsizenorm-{strand}.tsv")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}-{sample}-libsizenorm-{strand}.log"
    run:
        if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
        else:
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")

rule gzip_deeptools_matrix:
    input:
        tsv = "datavis/{annotation}/{annotation}-{sample}-libsizenorm-{strand}.tsv"
    output:
        "datavis/{annotation}/{annotation}-{sample}-libsizenorm-{strand}.tsv.gz"
    shell: """
        pigz -f {input.tsv}
        """

rule melt_matrix:
    input:
        matrix = "datavis/{annotation}/{annotation}-{sample}-libsizenorm-{strand}.tsv.gz"
    output:
        temp("datavis/{annotation}/{annotation}-{sample}-libsizenorm-{strand}-melted.tsv.gz")
    params:
        group = lambda wildcards : SAMPLES[wildcards.sample]["factor"] if wildcards.sample in SAMPLES else NEXUS_SAMPLES[wildcards.sample]["factor"],
        assay = lambda wildcards: "ChIP-exo" if wildcards.sample in SAMPLES else "ChIP-nexus",
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
    script:
        "scripts/melt_matrix.R"

rule cat_matrices:
    input:
        expand("datavis/{{annotation}}/{{annotation}}-{sample}-libsizenorm-{{strand}}-melted.tsv.gz", sample=list(SAMPLES.keys())+list(NEXUS_SAMPLES.keys()))
    output:
        "datavis/{annotation}/allsamples-{annotation}-libsizenorm-{strand}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{annotation}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule r_metagene:
    input:
        plus = "datavis/{annotation}/allsamples-{annotation}-libsizenorm-SENSE.tsv.gz",
        minus = "datavis/{annotation}/allsamples-{annotation}-libsizenorm-ANTISENSE.tsv.gz",
    params:
        trim_pct = lambda wildcards: config["annotations"][wildcards.annotation]["trim_pct"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
        downstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
        refptlabel= lambda wildcards: config["annotations"][wildcards.annotation]["refpointlabel"],
        factor = "TFIIB",
        ylabel = lambda wildcards: config["annotations"][wildcards.annotation]["ylabel"],
    output:
        "datavis/{annotation}/{annotation}-libsizenorm-exo-v-nexus-metagene.svg",
    script: "scripts/plotMetagene.R"
