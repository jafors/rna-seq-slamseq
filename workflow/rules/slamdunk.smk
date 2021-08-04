rule slamdunk_map:
    input:
        ref="resources/genome.fasta",
        files=get_reads,
    output:
        "results/map/{sample}_slamdunk_mapped.bam",
    params:
        extra=config["params"]["slamdunk"]["map"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        samplename=lambda wc: wc.sample,
        fastqname=lambda wc, input: os.path.basename(input.files[0]).replace(
            ".gz", ""
        ),
    log:
        "logs/slamdunk-map/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "slamdunk map -r {input.ref} -o {params.outdir} -sn {params.samplename} -t {threads} {params.extra} {input.files} > {log} 2>&1 "
        " && mv {params.outdir}/{params.fastqname}_slamdunk_mapped.bam {output}"


rule slamdunk_filter:
    input:
        bam="results/map/{sample}_slamdunk_mapped_sorted.bam",
        bai="results/map/{sample}_slamdunk_mapped_sorted.bam.bai",
        utr3="resources/utr3.bed",
    output:
        "results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
    params:
        extra=config["params"]["slamdunk"]["filter"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
    log:
        "logs/slamdunk-filter/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "slamdunk filter {input.bam} -o {params.outdir} -b {input.utr3} -t {threads} {params.extra} > {log} 2>&1"


rule slamdunk_snp:
    input:
        ref="resources/genome.fasta",
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
    output:
        "results/snp/{sample}_slamdunk_mapped_sorted_filtered_snp.vcf",
    params:
        extra=config["params"]["slamdunk"]["snp"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
    log:
        "logs/slamdunk-snp/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "slamdunk snp -r {input.ref} -o {params.outdir} -t {threads} {params.extra} {input.bam} > {log} 2>&1"


rule slamdunk_count:
    input:
        ref="resources/genome.fasta",
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
        snp="results/snp/{sample}_slamdunk_mapped_sorted_filtered_snp.vcf",
        utr3="resources/utr3.bed",
    output:
        "results/count/{sample}_slamdunk_mapped_sorted_filtered_tcount.tsv",
    params:
        extra=config["params"]["slamdunk"]["count"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"]),
    log:
        "logs/slamdunk-count/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "slamdunk count -r {input.ref} -o {params.outdir} -s {params.snp} -b {input.utr3} -t {threads} {params.extra} {input.bam} > {log} 2>&1"
