#!/usr/bin/env python

configfile: "config.yaml"

CSFASTA=config["exo_samples"]
EXO_SAMPLES = {k:v for e in [d for fq in [v["samples"] for k,v in CSFASTA.items()] for d in fq] for k,v in e.items()} #comprehension hell...just makes a dict where the top level keys are the samples
NEXUS_SAMPLES=config["nexus_samples"]
FACTORS = set(v["factor"] for (k,v) in EXO_SAMPLES.items())
INTERSECT_FACTORS = [x for x in FACTORS if x in set(v["factor"] for k,v in NEXUS_SAMPLES.items())]

localrules: all,
    make_barcode_fasta, cat_matrices

rule all:
    input:
        expand("csfasta/{sample}.csfastq", sample=EXO_SAMPLES),
        expand("csfasta/{sample}.csfastq", sample=list(EXO_SAMPLES.keys())),
        expand("peakcalling/macs/{factor}_peaks.narrowPeak", factor=FACTORS),
        #FastQC
        # expand("qual_ctrl/fastqc/{sample}/{sample}_fastqc.html", sample=list(SAMPLES.keys())+["unknown"]),
        expand("coverage/{norm}/{sample}-chipexo-{norm}-{strand}.{fmt}", norm=["counts","libsizenorm"], strand=["plus", "minus", "SENSE", "ANTISENSE", "protection"], sample=EXO_SAMPLES, fmt=["bedgraph", "bw"]),
        expand("correlations/{factor}-nexus-v-exo-window-{windowsize}-correlations.svg", factor=INTERSECT_FACTORS, windowsize=config["corr-windowsizes"]),
        expand("datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-{plottype}-bysample.svg", factor=INTERSECT_FACTORS, annotation=config["annotations"], plottype=["heatmap", "protection-metagene"])

rule make_barcode_fasta:
    output:
        "csfasta/{run}.fa"
    run:
        with open(str(output), "w") as out:
            for x in config["exo_samples"][wildcards.run]["samples"]:
                for k,v in x.items():
                    out.write(f'>{k}\n^{v["barcode"]}\n')

#minimum length 5 sanitizes output for bowtie
#output is uncompressed because I don't feel like learning how to open gzip files in perl right now...
#we do it in series rather than parallel because there was no straightforward way to specify the outputs otherwise
#note that 'unknown' reads are overwritten by each demultiplexing, and are not looked at further
rule trim_and_demultiplex:
    input:
        csfasta = [v["csfasta"] for k,v in CSFASTA.items()],
        qual = [v["qual"] for k,v in CSFASTA.items()],
        barcodes = expand("csfasta/{run}.fa", run=CSFASTA)
    params:
        trim_qual = config["cutadapt"]["trim_qual"]
    output:
        expand("csfasta/{sample}.csfastq", sample=list(EXO_SAMPLES.keys()))
    log: "logs/trim_and_demultiplex/trim_and_demultiplex.log"
    run:
        for run in list(CSFASTA.keys()):
            csfasta = CSFASTA[run]["csfasta"]
            qual = CSFASTA[run]["qual"]
            barcodes = "csfasta/" + run + ".fa"
            shell("""cutadapt -c -g file:{barcodes} -q {params.trim_qual} -m 5 -o csfasta/{{name}}.csfastq --length-tag 'length=' {csfasta} {qual};
                    mv csfasta/unknown.csfastq csfasta/unknown-{run}.csfastq""")

# rule fastqc:
#     input:
#         "csfasta/{sample}.csfastq"
#     output:
#         html = "qual_ctrl/fastqc/{sample}/{sample}_fastqc.html",
#         folder = "qual_ctrl/fastqc/{sample}/{sample}_fastqc.zip",
#     threads: config["threads"]
#     log: "logs/fastqc/fastqc-{sample}.log"
#     shell: """
#         (mkdir -p qual_ctrl/fastqc/{wildcards.sample}) &> {log}
#         (fastqc -o qual_ctrl/fastqc/{wildcards.sample} --noextract -t {threads} {input}) &>> {log}
#         """

# cutadapt outputs in fastq format, but bowtie throws an error when I try to align using the fastq
# (either "Reads file contained a pattern with more than 1024 quality values" using bowtie 1.2.1.1, or "Too few quality values for read: SRR..." with bowtie 1.1.2)
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

# #align with bowtie1 on colorspace index
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
    conda: "envs/bowtie.yaml"
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
    conda: "envs/bowtie.yaml"
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

rule macs2:
    input:
        bam = lambda wc: expand("alignment/{sample}.bam", sample={k for (k,v) in EXO_SAMPLES.items() if (v["factor"]==wc.factor)}),
        fasta = config["genome"]["fasta"]
    output:
        xls = "peakcalling/macs/{factor}_peaks.xls",
        peaks = "peakcalling/macs/{factor}_peaks.narrowPeak",
        summits = "peakcalling/macs/{factor}_summits.bed",
        script = "peakcalling/macs/{factor}_model.r",
        pdf = "peakcalling/macs/{factor}_model.pdf",
        treat_bg = "peakcalling/macs/{factor}_treat_pileup.bdg",
        cntrl_bg = "peakcalling/macs/{factor}_control_lambda.bdg"
    params:
        bw = config["macs2"]["bw"],
        slocal = config["macs2"]["slocal"],
        llocal = config["macs2"]["llocal"],
        qscore = config["macs2"]["fdr"]
    conda:
        "envs/macs2.yaml"
    log: "logs/macs2/macs2-{factor}.log"
    shell: """
        (macs2 callpeak -t {input.bam} -f BAM -g $(faidx {input.fasta} -i chromsizes | awk '{{sum += $2}} END {{print sum}}') --keep-dup all --bdg -n peakcalling/macs/{wildcards.factor} --SPMR --bw {params.bw} --slocal {params.slocal} --llocal {params.llocal} --call-summits -q {params.qscore}) &> {log}
        (Rscript peakcalling/macs/{wildcards.factor}_model.r) &>> {log}
        (sed -i -e 's/peakcalling\/macs\///g' peakcalling/macs/{wildcards.factor}_peaks.narrowPeak) &>> {log}
        (sed -i -e 's/peakcalling\/macs\///g' peakcalling/macs/{wildcards.factor}_summits.bed) &>> {log}
        """

rule get_protection:
    input:
        tsv = lambda wc: "peakcalling/macs/" + EXO_SAMPLES[wc.sample]["factor"] + "_peaks.xls",
        bam = "alignment/{sample}.bam"
    output:
        coverage = "coverage/counts/{sample}-chipexo-counts-protection.bedgraph"
    shell: """
        fragsize=$(grep -e "^# d = " {input.tsv} | cut -d ' ' -f4)
        genomeCoverageBed -bga -fs $fragsize -scale $(echo 1/$fragsize | bc -l) -ibam {input.bam} | LC_COLLATE=C sort -k1,1 -k2,2n > {output.coverage}
        """

rule normalize:
    input:
        counts = "coverage/counts/{sample}-chipexo-counts-{strand}.bedgraph",
        plmin = "coverage/counts/{sample}-chipexo-counts-combined.bedgraph",
    output:
        normalized = "coverage/libsizenorm/{sample}-chipexo-libsizenorm-{strand}.bedgraph",
    wildcard_constraints:
        strand="plus|minus|protection"
    log: "logs/normalize/normalize-{sample}-{strand}.log"
    shell: """
        (bash scripts/libsizenorm.sh {input.plmin} {input.counts} 1 > {output.normalized}) &> {log}
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

rule map_to_windows:
    input:
        bg = lambda wc: NEXUS_SAMPLES[wc.sample]["coverage"] + "SENSE.bedgraph" if wc.sample in NEXUS_SAMPLES else "coverage/libsizenorm/" + wc.sample + "-chipexo-libsizenorm-SENSE.bedgraph",
        fasta = config["genome"]["fasta"]
    output:
        exp = temp("coverage/{sample}-window-{windowsize}-coverage.bedgraph"),
    shell: """
        bedtools makewindows -g <(faidx {input.fasta} -i chromsizes | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}') -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output.exp}
        """

rule join_factor_window_counts:
    input:
        coverage = lambda wc: expand("coverage/{sample}-window-{windowsize}-coverage.bedgraph", sample=[k for k,v in NEXUS_SAMPLES.items() if v["factor"]==wc.factor] + [k for k,v in EXO_SAMPLES.items() if v["factor"]==wc.factor], windowsize=wc.windowsize)
    params:
        names = lambda wc: ["ChIP-nexus-"+str(idx+1) for idx,v in enumerate(NEXUS_SAMPLES.items()) if v[1]["factor"]==wc.factor] + ["ChIP-exo-"+str(idx+1) for idx,v in enumerate(EXO_SAMPLES.items()) if v[1]["factor"]==wc.factor]
    output:
        "correlations/union-bedgraph-{factor}-window-{windowsize}.tsv.gz"
    shell: """
        bedtools unionbedg -i {input.coverage} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output}
        """

rule plotcorrelations:
    input:
        "correlations/union-bedgraph-{factor}-window-{windowsize}.tsv.gz"
    output:
        "correlations/{factor}-nexus-v-exo-window-{windowsize}-correlations.svg"
    params:
        pcount = lambda wc: 0.01*int(wc.windowsize),
        samplelist = lambda wc: ["ChIP-nexus-"+str(idx+1) for idx,v in enumerate(NEXUS_SAMPLES.items()) if v[1]["factor"]==wc.factor] + ["ChIP-exo-"+str(idx+1) for idx,v in enumerate(EXO_SAMPLES.items()) if v[1]["factor"]==wc.factor]
    script:
        "scripts/plotcorr.R"

rule bedgraph_to_bigwig:
    input:
        bedgraph = "coverage/{norm}/{sample}-chipexo-{norm}-{strand}.bedgraph",
        fasta= config["genome"]["fasta"]
    output:
        "coverage/{norm}/{sample}-chipexo-{norm}-{strand}.bw",
    params:
        stranded = lambda wc: [] if wc.strand not in ["SENSE", "ANTISENSE"] else """| awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2; print $1"-minus", $2}}' | LC_COLLATE=C sort -k1,1"""
    log: "logs/bg_to_bw/bg_to_bw-{sample}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes {params.stranded}) {output}) &> {log}
        """

rule make_stranded_annotations:
    input:
        lambda wc : config["annotations"][wc.annotation]["path"]
    output:
        "datavis/{annotation}/{annotation}.bed"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (cut -f1-6 {input} | bash scripts/makeStrandedBed.sh > {output}) &> {log}
        """

rule compute_matrix:
    input:
        annotation = lambda wc: config["annotations"][wc.annotation]["path"] if wc.strand=="protection" else "datavis/{annotation}/{annotation}.bed".format(**wc),
        bw = lambda wc: "coverage/libsizenorm/" + wc.sample +"-chipexo-libsizenorm-"+ wc.strand + ".bw" if wc.sample in EXO_SAMPLES else NEXUS_SAMPLES[wc.sample]["coverage"] + wc.strand + ".bw"
    output:
        dtfile = temp("datavis/{annotation}/libsizenorm/{annotation}_{sample}_{factor}-libsizenorm-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/libsizenorm/{annotation}_{sample}_{factor}-libsizenorm-{strand}.tsv"),
        melted = "datavis/{annotation}/libsizenorm/{annotation}_{sample}_{factor}-libsizenorm-{strand}-melted.tsv.gz",
    params:
        group = lambda wc : "ChIP-exo" if wc.sample in EXO_SAMPLES else "ChIP-nexus",
        upstream = lambda wc: config["annotations"][wc.annotation]["upstream"] + config["annotations"][wc.annotation]["binsize"],
        dnstream = lambda wc: config["annotations"][wc.annotation]["dnstream"] + config["annotations"][wc.annotation]["binsize"],
        binsize = lambda wc: config["annotations"][wc.annotation]["binsize"],
        sort = lambda wc: config["annotations"][wc.annotation]["sort"],
        sortusing = lambda wc: config["annotations"][wc.annotation]["sortby"],
        binstat = lambda wc: config["annotations"][wc.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/compute_matrix/compute_matrix-{annotation}_{sample}_{factor}-libsizenorm-{strand}.log"
    run:
        if config["annotations"][wildcards.annotation]["type"]=="absolute":
            refpoint = config["annotations"][wildcards.annotation]["refpoint"]
            if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
                shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
            else:
                shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            refpoint = "TSS"
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {refpoint} --group {params.group} -s {wildcards.sample} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wc: expand("datavis/{annotation}/libsizenorm/{annotation}_{sample}_{factor}-libsizenorm-{strand}-melted.tsv.gz", sample=[k for k,v in EXO_SAMPLES.items() if v["factor"]==wc.factor] + [k for k,v in NEXUS_SAMPLES.items() if v["factor"]==wc.factor], annotation=wc.annotation, factor=wc.factor, strand=wc.strand)
    output:
        "datavis/{annotation}/libsizenorm/allsamples-{annotation}-{factor}-libsizenorm-{strand}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{annotation}_{factor}-libsizenorm-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_heatmaps:
    input:
        matrix = "datavis/{annotation}/libsizenorm/allsamples-{annotation}-{factor}-libsizenorm-protection.tsv.gz"
    output:
        heatmap_sample = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-heatmap-bysample.svg",
        heatmap_group = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-heatmap-bygroup.svg",
    params:
        samplelist = lambda wc: [k for k,v in EXO_SAMPLES.items() if v["factor"]==wc.factor] + [k for k,v in NEXUS_SAMPLES.items() if v["factor"]==wc.factor],
        mtype = lambda wc : config["annotations"][wc.annotation]["type"],
        upstream = lambda wc : config["annotations"][wc.annotation]["upstream"],
        dnstream = lambda wc : config["annotations"][wc.annotation]["dnstream"],
        pct_cutoff = lambda wc : config["annotations"][wc.annotation]["pct_cutoff"],
        cluster = lambda wc : config["annotations"][wc.annotation]["cluster"],
        nclust = lambda wc: config["annotations"][wc.annotation]["nclusters"],
        heatmap_cmap = lambda wc : config["annotations"][wc.annotation]["heatmap_colormap"],
        refpointlabel = lambda wc : config["annotations"][wc.annotation]["refpointlabel"],
        ylabel = lambda wc : config["annotations"][wc.annotation]["ylabel"]
    run:
        if config["annotations"][wildcards.annotation]["type"]=="scaled":
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            endlabel = config["annotations"][wildcards.annotation]["three_prime_label"]
        else:
            scaled_length=0
            endlabel = "HAIL SATAN!"
        shell("""Rscript scripts/plot_nexus_heatmaps.R -i {input.matrix} -s {params.samplelist} -t {params.mtype} -u {params.upstream} -d {params.dnstream} -c {params.pct_cutoff} -z {params.cluster} -k {params.nclust} -r {params.refpointlabel} -f {wildcards.factor} -l {scaled_length} -e {endlabel} -y {params.ylabel} -m {params.heatmap_cmap} -o {output.heatmap_sample} -p {output.heatmap_group}""")

rule plot_metagenes:
    input:
        plus = "datavis/{annotation}/libsizenorm/allsamples-{annotation}-{factor}-libsizenorm-SENSE.tsv.gz",
        minus = "datavis/{annotation}/libsizenorm/allsamples-{annotation}-{factor}-libsizenorm-ANTISENSE.tsv.gz",
        protection = "datavis/{annotation}/libsizenorm/allsamples-{annotation}-{factor}-libsizenorm-protection.tsv.gz"
    output:
        smeta_group = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-stranded-metagene-bygroup.svg",
        smeta_sample = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-stranded-metagene-bysample.svg",
        pmeta_group = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-protection-metagene-bygroup.svg",
        pmeta_sample = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-protection-metagene-bysample.svg",
        pmeta_goverlay = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-protection-metagene-overlay-group.svg",
        pmeta_soverlay = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-protection-metagene-overlay-sample.svg",
        pmeta_soverlay_bygroup = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-protection-metagene-overlay-sample-bygroup.svg",
        meta_heatmap_group = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-protection-metaheatmap-bygroup.svg",
        meta_heatmap_sample = "datavis/{annotation}/libsizenorm/{factor}-{annotation}-libsizenorm-nexus-v-exo-protection-metaheatmap-bysample.svg",
    params:
        samplelist = lambda wc: [k for k,v in EXO_SAMPLES.items() if v["factor"]==wc.factor] + [k for k,v in NEXUS_SAMPLES.items() if v["factor"]==wc.factor],
        mtype = lambda wc : config["annotations"][wc.annotation]["type"],
        upstream = lambda wc : config["annotations"][wc.annotation]["upstream"],
        dnstream = lambda wc : config["annotations"][wc.annotation]["dnstream"],
        trim_pct = lambda wc : config["annotations"][wc.annotation]["trim_pct"],
        refpointlabel = lambda wc : config["annotations"][wc.annotation]["refpointlabel"],
        ylabel = lambda wc : config["annotations"][wc.annotation]["ylabel"]
    run:
        if config["annotations"][wildcards.annotation]["type"]=="scaled":
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            endlabel = config["annotations"][wildcards.annotation]["three_prime_label"]
        else:
            scaled_length=0
            endlabel = "HAIL SATAN!"
        shell("""Rscript scripts/plot_nexus_metagenes.R --inplus {input.plus} --inminus {input.minus} --inprotection {input.protection} -s {params.samplelist} -t {params.mtype} -u {params.upstream} -d {params.dnstream} -p {params.trim_pct} -r {params.refpointlabel} -f {wildcards.factor} -l {scaled_length} -e {endlabel} -y {params.ylabel} --out1 {output.smeta_group} --out2 {output.smeta_sample} --out3 {output.pmeta_group} --out4 {output.pmeta_sample} --out5 {output.pmeta_goverlay} --out6 {output.pmeta_soverlay} --out7 {output.pmeta_soverlay_bygroup} --out8 {output.meta_heatmap_group} --out9 {output.meta_heatmap_sample}""")

