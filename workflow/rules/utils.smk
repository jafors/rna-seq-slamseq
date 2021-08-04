rule sambamba:
    input:
        "{prefix}.bam",
    output:
        "{prefix}_sorted.bam",
        "{prefix}_sorted.bam.bai",
    log:
        "logs/sambamba-sort/{prefix}.log",
    wrapper:
        "0.77.0/bio/sambamba/sort"
