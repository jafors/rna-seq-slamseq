rule dedup:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
    output:
        "results/filter/{sample}_slamdunk_mapped_sorted_filtered_dedup.bam",
    params:
        extra=config["params"]["alleyoop"]["dedup"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
    log:
        "logs/alleyoop-dedup/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 5
    shell:
        "alleyoop dedup -o {params.outdir} -t {threads} {params.extra} {input.bam} > {log} 2>&1"


rule collapse:
    input:
        tcount="results/count/{sample}_slamdunk_mapped_sorted_filtered_tcount.tsv",
    output:
        "results/count/{sample}_slamdunk_mapped_sorted_filtered_tcount_collapsed.tsv",
    params:
        extra=config["params"]["alleyoop"]["collapse"],
        outdir=lambda wc, input, output: os.path.dirname(output),
    log:
        "logs/alleyoop-collapse/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 5
    shell:
        "alleyoop collapse -o {params.outdir} -t {threads} {params.extra} {input.tcount} > {log} 2>&1"


rule rates:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
        ref="resources/genome.fasta",
    output:
        multiext(
            "results/alleyoop/{sample}_slamdunk_mapped_sorted_filtered_overallrates",
            ".csv",
            ".pdf",
        ),
    params:
        extra=config["params"]["alleyoop"]["rates"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
    log:
        "logs/alleyoop-rates/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop rates -r {input.ref} -o {params.outdir} -t {threads} {params.extra} {input.bam} > {log} 2>&1"


rule tccontext:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
        ref="resources/genome.fasta",
    output:
        multiext(
            "results/alleyoop/{sample}_slamdunk_mapped_sorted_filtered_tccontext",
            ".csv",
            ".pdf",
        ),
    params:
        extra=config["params"]["alleyoop"]["tccontext"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
    log:
        "logs/alleyoop-tccontext/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop tccontext -r {input.ref} -o {params.outdir} -t {threads} {params.extra} {input.bam} > {log} 2>&1"


rule utrrates:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
        ref="resources/genome.fasta",
        utr3="resources/utr3.bed",
    output:
        multiext(
            "results/alleyoop/{sample}_slamdunk_mapped_sorted_filtered_mutationrates_utr",
            ".csv",
            ".pdf",
        ),
    params:
        extra=config["params"]["alleyoop"]["utrrates"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
    log:
        "logs/alleyoop-utrrates/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop utrrates -r {input.ref} -b {input.utr3} -o {params.outdir} -t {threads} {params.extra} {input.bam} > {log} 2>&1"


rule snpeval:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
        ref="resources/genome.fasta",
        utr3="resources/utr3.bed",
        snp="results/snp/{sample}_slamdunk_mapped_sorted_filtered_snp.vcf",
    output:
        multiext(
            "results/alleyoop/{sample}_slamdunk_mapped_sorted_filtered_SNPeval",
            ".csv",
            ".pdf",
        ),
    params:
        extra=config["params"]["alleyoop"]["snpeval"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"]),
    log:
        "logs/alleyoop-snpeval/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop snpeval -r {input.ref} -b {input.utr3} -s {params.snp} -o {params.outdir} -t {threads} {params.extra} {input.bam} > {log} 2>&1"


rule summary:
    input:
        tcount=expand(
            "results/count/{sample}_slamdunk_mapped_sorted_filtered_tcount.tsv",
            sample=samples["sample_name"],
        ),
        bam=expand(
            "results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
            sample=samples["sample_name"],
        ),
    output:
        multiext("results/alleyoop/summary", ".txt", "_PCA.txt", "_PCA.pdf"),
    log:
        "logs/alleyoop-summary.log",
    params:
        countdir=lambda wc, input: os.path.dirname(input["tcount"][0]),
        outfile=lambda wc, input, output: output[0],
    conda:
        "../envs/slamdunk.yaml"
    shell:
        "alleyoop summary -o {params.outfile} -t {params.countdir} {input.bam} > {log} 2>&1"


rule tcperreadpos:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
        ref="resources/genome.fasta",
        snp="results/snp/{sample}_slamdunk_mapped_sorted_filtered_snp.vcf",
    output:
        multiext(
            "results/alleyoop/{sample}_slamdunk_mapped_sorted_filtered_tcperreadpos",
            ".csv",
            ".pdf",
        ),
    params:
        extra=config["params"]["alleyoop"]["tcperreadpos"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"]),
    log:
        "logs/alleyoop-tcperreadpos/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop tcperreadpos -r {input.ref} -o {params.outdir} -s {params.snp} -t {threads} {params.extra} {input.bam} > {log} 2>&1"


rule tcperutrpos:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
        ref="resources/genome.fasta",
        utr3="resources/utr3.bed",
        snp="results/snp/{sample}_slamdunk_mapped_sorted_filtered_snp.vcf",
    output:
        multiext(
            "results/alleyoop/{sample}_slamdunk_mapped_sorted_filtered_tcperutr",
            ".csv",
            ".pdf",
        ),
    params:
        extra=config["params"]["alleyoop"]["tcperutrpos"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"]),
    log:
        "logs/alleyoop-tcperutrpos/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    threads: 10
    shell:
        "alleyoop tcperutrpos -r {input.ref} -o {params.outdir} -b {input.utr3} -s {params.snp} -t {threads} {params.extra} {input.bam} > {log} 2>&1"


rule dump:
    input:
        bam="results/filter/{sample}_slamdunk_mapped_sorted_filtered.bam",
        ref="resources/genome.fasta",
        snp="results/snp/{sample}_slamdunk_mapped_sorted_filtered_snp.vcf",
    output:
        "results/alleyoop/{sample}_slamdunk_mapped_sorted_filtered_readinfo.sdunk",
    params:
        extra=config["params"]["alleyoop"]["dump"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"]),
    log:
        "logs/alleyoop-dump/{sample}.log",
    conda:
        "../envs/slamdunk.yaml"
    shell:
        "alleyoop dump -r {input.ref} -o {params.outdir} -s {params.snp} {params.extra} {input.bam} > {log} 2>&1"
