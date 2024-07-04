rule deconvolutexengsort:
   input:
      f1 =  get_fq1,
      f2 =  get_fq2
   output:
      graftf1 =  "results/xengsort/{sample}-graft.1.fq.gz",
      graftf2 =  "results/xengsort/{sample}-graft.2.fq.gz"
    #neitherf1 =  temp("results/xengsort/{sample}-neither.1.fq"),
    #neitherf2 =  temp("results/xengsort/{sample}-neither.2.fq"),
    #bothf1 =  temp("results/xengsort/{sample}-both.1.fq"),
    #bothf2 =  temp("results/xengsort/{sample}-both.2.fq"),
    #ambiguousf1 =  temp("results/xengsort/{sample}-ambiguous.1.fq"),
    #ambiguousf2 =  temp("results/xengsort/{sample}-ambiguous.2.fq"),
    #hostf1 =  temp("results/xengsort/{sample}-host.1.fq"),
    #hostf2 =  temp("results/xengsort/{sample}-host.2.fq")
   params:
      xengsortidx = config["ref"]["xengsortidx"],
      xengsortcontainer = config['env']['xengsort'],
      sampleid = "{sample}"
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
        alignedcoord="results/star/pe/{sample}/Aligned.sortedByCoord.out.bam",
        alignedtranscriptome="results/star/pe/{sample}/Aligned.toTranscriptome.out.bam",
        counts="results/star/pe/{sample}/ReadsPerGene.out.tab",
    threads: 16
    params:
        stargtf=config['star']['gtf'],
        starparams=config['star']['params'],
        stargenome=config["star"]["star-genome"],
        outprefix="results/star/pe/{sample}/",
    shell:
        """
        module load STAR/2.7.3a
    
        STAR \
        --quantMode GeneCounts TranscriptomeSAM \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterIntronMotifs RemoveNoncanonical \
        --chimSegmentMin 10 \
        --chimOutType SeparateSAMold \
        --outSAMunmapped Within \
        --sjdbGTFfile {params.stargtf} {params.starparams} \
        --runThreadN {threads} \
        --genomeDir {params.stargenome} \
        --readFilesIn {input.fq1} {input.fq2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.outprefix} \
        --outStd Log
        """

rule align_se:
    input:
        fq1="results/xengsort/{sample}-graft.1.fq.gz",
    output:
      alignedcoord="results/star/se/{sample}/Aligned.sortedByCoord.out.bam",
      alignedtranscriptome="results/star/se/{sample}/Aligned.toTranscriptome.out.bam",
      counts="results/star/se/{sample}/ReadsPerGene.out.tab",
    params:
      stargtf=config['star']['gtf'],
      starparams=config['star']['params'],
      stargenome=config["star"]["star-genome"],
      outprefix="results/star/se/{sample}/",
    threads: 16
    shell:
      """
      module load STAR/2.7.3a
      
      STAR \
      --quantMode GeneCounts TranscriptomeSAM \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterIntronMotifs RemoveNoncanonical \
      --chimSegmentMin 10 \
      --chimOutType SeparateSAMold \
      --outSAMunmapped Within \
      --sjdbGTFfile {params.stargtf} {params.starparams} \
      --runThreadN {threads} \
      --genomeDir {params.stargenome} \
      --readFilesIn {input.fq1} \
      --readFilesCommand zcat \
      --outFileNamePrefix {params.outprefix} \
      --outStd Log
      """

rule index_coord:
    input:
      coord="results/star/{strand}/{sample}/Aligned.sortedByCoord.out.bam",
    output:
      coord="results/star/{strand}/{sample}/Aligned.sortedByCoord.out.bam.bai",
    params:
    shell:
      """
      module load samtools/1.17
      
      samtools index {input.coord} -o {output.coord}
      """

rule index_coord_transcriptome:
    input:
      transcriptome="results/star/{strand}/{sample}/Aligned.toTranscriptome.out.bam",
    output:
      transcriptome="results/star/{strand}/{sample}/Aligned.toTranscriptome.out.bam.bai",
    params:
    shell:
      """
      module load samtools/1.17
    
      samtools index {input.transcriptome} -o {output.transcriptome}
      """