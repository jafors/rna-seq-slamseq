rule get_genome:
    output:
        "resources/genome.fasta"
    log:
        "logs/get-genome.log"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    cache: True
    wrapper:
        "0.59.2/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf"
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="" # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    cache: True
    log:
        "logs/get_annotation.log"
    wrapper:
        "0.59.2/bio/reference/ensembl-annotation"

rule gtf_to_bed:
    input:
        "resources/genome.gtf"
    output:
        "resources/utr3.bed"
    log:
        "logs/gtf2bed.log"
    conda:
        "bedops.yaml"
    shell:
        "grep 'three_prime_utr' {input} | gtf2bed - > {output}"


rule genome_faidx:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.fasta.fai"
    log:
        "logs/genome-faidx.log"
    cache: True
    wrapper:
        "0.59.2/bio/samtools/faidx"