rule dedup:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam"
    output:
        "results/filter/{sample}_slamdunk_mapped_filtered_dedup.bam"
    params:
        extra=config["params"]["alleyoop"]["dedup"],
        outdir=lambda wc, input, output: os.path.dirname(output[0])
    log:
        "logs/alleyoop-dedup/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 5
    shell:
        "alleyoop dedup -o {params.outdir} -t {threads} {params.extra} {input.bam}"

rule collapse:
    input:
        tcount="results/count/{sample}_slamdunk_mapped_filtered_tcount.tsv"
    output:
        "results/count/{sample}_slamdunk_mapped_filtered_tcount_collapsed.tsv"
    params:
        extra=config["params"]["alleyoop"]["collapse"],
        outdir=lambda wc, input, output: os.path.dirname(output[0])
    log:
        "logs/alleyoop-collapse/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 5
    shell:
        "alleyoop collapse -o {params.outdir} -t {threads} {params.extra} {input.tcount}"

rule rates:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam",
        ref="resources/genome.fasta"
    output:
        multiext("results/rates/{sample}_slamdunk_mapped_filtered_overallrates", ".csv", ".pdf")
    params:
        extra=config["params"]["alleyoop"]["rates"],
        outdir=lambda wc, input, output: os.path.dirname(output[0])
    log:
        "logs/alleyoop-rates/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop rates -r {input.ref} -o {params.outdir} -t {threads} {params.extra} {input.bam}"

rule tccontext:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam",
        ref="resources/genome.fasta"
    output:
        multiext("results/tccontext/{sample}_slamdunk_mapped_filtered_tccontext", ".csv", ".pdf")
    params:
        extra=config["params"]["alleyoop"]["tccontext"],
        outdir=lambda wc, input, output: os.path.dirname(output[0])
    log:
        "logs/alleyoop-tccontext/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop rates -r {input.ref} -o {params.outdir} -t {threads} {params.extra} {input.bam}"

rule utrrates:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam",
        ref="resources/genome.fasta",
        utr3="resources/utr3.bed"
    output:
        multiext("results/utrrates/{sample}_slamdunk_mapped_filtered_mutationrates_utr", ".csv", ".pdf")
    params:
        extra=config["params"]["alleyoop"]["utrrates"],
        outdir=lambda wc, input, output: os.path.dirname(output[0])
    log:
        "logs/alleyoop-utrrates/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop rates -r {input.ref} -b {input.utr3} -o {params.outdir} -t {threads} {params.extra} {input.bam}"

rule snpeval:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam",
        ref="resources/genome.fasta",
        utr3="resources/utr3.bed",
        snp="results/snp/{sample}_slamdunk_mapped_filtered.vcf"
    output:
        multiext("results/snpeval/{sample}_slamdunk_mapped_filtered_SNPeval", ".csv", ".pdf")
    params:
        extra=config["params"]["alleyoop"]["snpeval"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"])
    log:
        "logs/alleyoop-snpeval/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop rates -r {input.ref} -b {input.utr3} -s {params.snp} -o {params.outdir} -t {threads} {params.extra} {input.bam}"

rule summary:
    input:
        tcount=expand("results/count/{sample}_slamdunk_mapped_filtered_tcount.tsv", sample=samples["sample_name"]),
        bam=expand("results/filter/{sample}_slamdunk_mapped_filtered.bam", sample=samples["sample_name"])
    output:
        multiext("results/summary", ".txt", "_PCA.txt", "_PCA.pdf")
    log:
        "logs/alleyoop-summary.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 5
    shell:
        "alleyoop summary -o {output[0]} -t {threads} -t {input.tcount} {input.bam}"

rule tcperreadpos:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam",
        ref="resources/genome.fasta",
        snp="results/snp/{sample}_slamdunk_mapped_filtered.vcf"
    output:
        multiext("results/{sample}_slamdunk_mapped_filtered_tcperreadpos", ".csv", ".pdf")
    params:
        extra=config["params"]["alleyoop"]["tcperreadpos"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"])
    log:
        "logs/alleyoop-tcperreadpos/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop tcperreadpos -r {input.ref} -o {params.outdir} -s {params.snp} -t {threads} {params.extra} {input.bam}"
    
rule tcperutrpos:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam",
        ref="resources/genome.fasta",
        utr3="resources/utr3.bed",
        snp="results/snp/{sample}_slamdunk_mapped_filtered.vcf"
    output:
        multiext("results/{sample}_slamdunk_mapped_filtered_tcperutr", ".csv", ".pdf")
    params:
        extra=config["params"]["alleyoop"]["tcperutrpos"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"])
    log:
        "logs/alleyoop-tcperutrpos/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop tcperutrpos -r {input.ref} -o {params.outdir} -b {input.utr3} -s {params.snp} -t {threads} {params.extra} {input.bam}"

rule dump:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam",
        ref="resources/genome.fasta",
        snp="results/snp/{sample}_slamdunk_mapped_filtered.vcf"
    output:
        "results/{sample}_slamdunk_mapped_filtered_readinfo.sdunk"
    params:
        extra=config["params"]["alleyoop"]["dump"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"])
    log:
        "logs/alleyoop-dump/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    shell:
        "alleyoop dump -r {input.ref} -o {params.outdir} -s {params.snp} {params.extra} {input.bam}"
