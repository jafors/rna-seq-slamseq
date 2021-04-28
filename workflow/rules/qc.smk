rule fastqc:
    input:
        get_reads
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 1
    wrapper:
        "0.74.0/bio/fastqc"

rule multiqc:
    input:
        get_multiqc_input,
        get_fastqc,
        "results/summary.txt",
    output:
        "results/qc/multiqc/multiqc.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc results -o results/qc/multiqc"
