rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "0.59.2/bio/sra-tools/fasterq-dump"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/catadapt/{sample}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="results/trimmed/{sample}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}_R2.fastq.gz",
        qc="results/trimmed/{sample}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}.log",
    params:
        others=config["params"]["cutadapt-pe"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "0.59.2/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="results/trimmed/{sample}_single.fastq.gz",
        qc="results/trimmed/{sample}_single.qc.txt",
    log:
        "logs/cutadapt/{sample}.log",
    params:
        others=config["params"]["cutadapt-se"],
        adapters_r1=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "0.59.2/bio/cutadapt/se"


rule merge_fastqs:
    input:
        get_fastqs,
    output:
        "results/merged/{sample}_{read}.fastq.gz",
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    shell:
        "cat {input} > {output} 2> {log}"