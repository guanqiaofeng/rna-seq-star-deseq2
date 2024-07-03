rule deconvolutexengsort:
   input:
      f1 =  get_map_reads_input_R1,
      f2 =  get_map_reads_input_R2
   output:
      graftf1 =  "results/xengsort/{sample}-graft.1.fq.gz",
      graftf2 =  "results/xengsort/{sample}-graft.2.fq.gz",
    #neitherf1 =  temp("results/xengsort/{sample}-neither.1.fq"),
    #neitherf2 =  temp("results/xengsort/{sample}-neither.2.fq"),
    #bothf1 =  temp("results/xengsort/{sample}-both.1.fq"),
    #bothf2 =  temp("results/xengsort/{sample}-both.2.fq"),
    #ambiguousf1 =  temp("results/xengsort/{sample}-ambiguous.1.fq"),
    #ambiguousf2 =  temp("results/xengsort/{sample}-ambiguous.2.fq"),
    #hostf1 =  temp("results/xengsort/{sample}-host.1.fq"),
    #hostf2 =  temp("results/xengsort/{sample}-host.2.fq")
   params:
      xengsortidx=config["ref"]["xengsortidx"],
      xengsortcontainer=config['env']['xengsort'],
      sampleid="{sample}"
   threads: 4
   shell:
      """
      module load apptainer/1.0.2
      module load pigz/2.6 
      
      mkdir -p results/xengsort
      mkdir tmp
      zcat {input.f1} > tmp/{params.sampleid}_R1.fastq
      zcat {input.f2} > tmp/{params.sampleid}_R2.fastq
      
      apptainer run {params.xengsortcontainer} \
      xengsort classify \
      --index {params.xengsortidx} \
      --fastq tmp/{params.sampleid}_R1.fastq \
      --pairs tmp/{params.sampleid}_R2.fastq \
      --prefix results/xengsort/{params.sampleid} \
      --compression none \
      -T {threads} \
      --progress \
      --filter
      
      rm tmp/{params.sampleid}_R1.fastq
      rm tmp/{params.sampleid}_R2.fastq
      pigz -{threads} results/xengsort/{params.sampleid}-graft.1.fq
      pigz -{threads} results/xengsort/{params.sampleid}-graft.2.fq
      """

rule align_pe:
    input:
        fq1="results/xengsort/{sample}-graft.1.fq.gz",
        fq2="results/xengsort/{sample}-graft.2.fq.gz",
    output:
        "results/star/pe/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/pe/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        "results/star/pe/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star-pe/{sample}-{unit}.log",
    params:
        index=config["star"]["star-genome"],
        extra="--quantMode GeneCounts TranscriptomeSAM "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterIntronMotifs RemoveNoncanonical "
        "--chimSegmentMin 10 "
        "--chimOutType SeparateSAMold "
        "--outSAMunmapped Within "
        "--sjdbGTFfile {} {}".format(
            config["star"]["gtf"], config["star"]["params"]
        ),
    threads: 24
    wrapper:
        "v0.75.0/bio/star/align"

rule align_se:
    input:
        fq1="results/xengsort/{sample}-graft.1.fq.gz",
    output:
        "results/star/se/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/se/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        "results/star/se/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star-se/{sample}-{unit}.log",
    params:
        index=config["star"]["star-genome"],
        extra="--quantMode GeneCounts TranscriptomeSAM "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterIntronMotifs RemoveNoncanonical "
        "--chimSegmentMin 10 "
        "--chimOutType SeparateSAMold "
        "--outSAMunmapped Within "
        "--sjdbGTFfile {} {}".format(
            config["star"]["gtf"], config["star"]["params"]
        ),
    threads: 24
    wrapper:
        "v0.75.0/bio/star/align"

rule index_coord:
  input:
    get_star_bam,
  output:
    "results/star/{ends}/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai",
  log:
    "logs/samtools/index/{sample}-{unit}.{ends}.sortedByCoord.log"
  wrapper:
    "v0.75.0/bio/samtools/index"
