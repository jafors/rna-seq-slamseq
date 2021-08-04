from snakemake.utils import validate
import pandas as pd
import os
import itertools

##### load config and sample sheets #####


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)


def get_multiqc_input(wildcards):
    # data = multiext(, "_tcperutr.csv", "tcperreadpos.csv"),
    data = list(
        itertools.chain.from_iterable(
            [
                multiext(
                    x,
                    "_tcperutr.csv",
                    "_tcperreadpos.csv",
                    "_overallrates.csv",
                    "_tccontext.csv",
                    "_mutationrates_utr.csv",
                    "_SNPeval.csv",
                )
                for x in expand(
                    "results/alleyoop/{sample}_slamdunk_mapped_sorted_filtered",
                    sample=samples["sample_name"],
                )
            ]
        )
    )
    return data


def get_fastqc(wildcards):
    return expand(
        "results/qc/fastqc/{sample}_fastqc.zip", sample=samples["sample_name"]
    )


def get_reads(wildcards):
    sample_units = units.loc[wildcards.sample]
    return sample_units["fq1"]
