def get_input_files(wildcards):
    sample_units = units.loc[wildcards.sample]
    return sample_units["fq1"]

rule slamdunk_map:
    input:
        ref="resources/genome.fasta",
        files=get_input_files
    output:
        "results/map/{sample}_slamdunk_mapped.bam"
    params:
        extra=config["params"]["slamdunk"]["map"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        samplename= lambda wc: wc.sample
    log:
        "logs/slamdunk-map/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads:
        10
    shell:
        "slamdunk map -r {input.ref} -o {params.outdir} -sn {params.samplename} -t {threads} {params.extra} {input.files}"

rule slamdunk_filter:
    input:
        ref="resources/genome.fasta",
        bam="results/map/{sample}_slamdunk_mapped.bam",
        utr3="resources/utr3.bed"
    output:
        "results/filter/{sample}_slamdunk_mapped_filtered.bam"
    params:
        extra=config["params"]["slamdunk"]["filter"],
        outdir=lambda wc, input, output: os.path.dirname(output[0])
    log:
        "logs/slamdunk-filter/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads:
        10
    shell:
        "slamdunk filter -r {input.ref} -o {params.outdir} -b {input.utr3} -t {threads} {params.extra} {input.bam}"

rule slamdunk_snp:
    input:
        ref="resources/genome.fasta",
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam"
    output:
        "results/snp/{sample}_slamdunk_mapped_filtered.vcf"
    params:
        extra=config["params"]["slamdunk"]["snp"],
        outdir=lambda wc, input, output: os.path.dirname(output[0])
    log:
        "logs/slamdunk-snp/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads:
        10
    shell:
        "slamdunk snp -r {input.ref} -o {output} -t {threads} {params.extra} {input.bam}"

rule slamdunk_count:
    input:
        ref="resources/genome.fasta",
        bam="results/filter/{sample}_slamdunk_mapped_filtered.bam",
        snp="results/snp/{sample}_slamdunk_mapped_filtered.vcf",
        utr3="resources/utr3.bed"
    output:
        "results/count/{sample}_slamdunk_mapped_filtered_tcount.tsv"
    params:
        extra=config["params"]["slamdunk"]["count"],
        outdir=lambda wc, input, output: os.path.dirname(output[0]),
        snp=lambda wc, input: os.path.dirname(input["snp"])
    log:
        "logs/slamdunk-count/{sample}.log"
    conda:
        "../envs/slamdunk.yaml"
    threads:
        10
    shell:
        "slamdunk count -r {input.ref} -o {params.outdir} -s {params.snp} -b {input.utr3} -t {threads} {params.extra} {input.bam}"