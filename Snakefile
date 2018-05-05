#from os import path
import glob

#Logan = 'AVEAY010B'.split()
#Buri = 'AVEAY001A'.split()
#Buri_calf  = 'AVEAY001B AVEAY002A AVEAY002B AVEAY003A AVEAY003B AVEAY004A'.split()
#Buri_calf_dam = 'AVEAY004B AVEAY005A AVEAY005B AVEAY006A AVEAY006B AVEAY007A'.split()
#HH_calf = 'AVEAY007B AVEAY008A AVEAY008B'.split()
#HH_calf_dam = 'AVEAY009A AVEAY009B AVEAY010A'.split()
#SAMPLES = Logan + Buri + Buri_calf + Buri_calf_dam + HH_calf + HH_calf_dam

fastq_pattern = "data/fastq/*/*_001.fastq.gz"
fastq_files = glob.glob(fastq_pattern)
sample_ids = [fq[11:-9] for fq in fastq_files]
fragment_ids = [fq[11:-16] for fq in fastq_files]

rule all:
    input:
#        expand('output/{sample}.txt', sample=SAMPLES)#, R=['R1', 'R2'])
        expand("qc/fastq/{sample}_fastqc.html", sample=sample_ids), expand("qc/fastq/{sample}_fastqc.zip", sample=sample_ids),
        "multiqc/fastq/multiqc_report.html",
        expand('data/trimmed/{fragment}_R1_001.fastq.gz', fragment=fragment_ids), expand('data/trimmed/{fragment}_R2_001.fastq.gz', fragment=fragment_ids),
        expand("qc/trimmed/{sample}_fastqc.html", sample=sample_ids), expand("qc/trimmed/{sample}_fastqc.zip", sample=sample_ids),
        "multiqc/trimmed/multiqc_report.html",
        "refGenome/ARS-UCD1.2_chr.fa",
        expand("refGenome/BwaIndex/genome.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
        expand("data/mapped_reads/{fragment}.bam", fragment=fragment_ids),

rule fastqc_pre:
    input:
        "data/fastq/{sample}.fastq.gz"
    output:
        html="qc/fastq/{sample}_fastqc.html",
        zip="qc/fastq/{sample}_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    threads:1
    shell:
        '''
        module load FastQC/0.11.5
        fastqc -t {threads} -f fastq -noextract -o qc/fastq/$(dirname {wildcards.sample}) {input}
        '''

rule multiQC_pre:
    input:
        html=expand("qc/fastq/{sample}_fastqc.html", sample=sample_ids)
    output:
        "multiqc/fastq/multiqc_report.html",
        "multiqc/fastq/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        "/opt/software/multiQC/1.0--singularity/bin/multiqc.img -z -o multiqc/fastq qc/fastq" 

rule trimmomatic_pe:
    input:
        r1="data/fastq/{fragment}_R1_001.fastq.gz",
        r2="data/fastq/{fragment}_R2_001.fastq.gz"
    output:
        r1="data/trimmed/{fragment}_R1_001.fastq.gz",
        r2="data/trimmed/{fragment}_R2_001.fastq.gz",
        sr1="data/trimmed/{fragment}_R1.unpaired_001.fastq.gz",
        sr2="data/trimmed/{fragment}_R2.unpaired_001.fastq.gz"
#    log:
#        "logs/trimmomatic/{fragment}.log"
#    params:
#        # list of trimmers (see manual)
#        trimmer=["TRAILING:3"],
#        # optional parameters
#        extra=""
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 4
    threads:2
    shell:
        '''
        module load Trimmomatic/0.33
        java -jar $TRIM/trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} {output.r1} {output.sr1} {output.r2} {output.sr2} ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:20
        '''

rule fastqc_post:
    input:
        "data/trimmed/{sample}.fastq.gz"
    output:
        html="qc/trimmed/{sample}_fastqc.html",
        zip="qc/trimmed/{sample}_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    threads:1
    shell:
        '''
        module load FastQC/0.11.5
        fastqc -t {threads} -f fastq -noextract -o qc/trimmed/$(dirname {wildcards.sample}) {input}
        '''

rule multiQC_post:
    input:
        html=expand("qc/trimmed/{sample}_fastqc.html", sample=sample_ids)
    output:
        "multiqc/trimmed/multiqc_report.html",
        "multiqc/trimmed/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        "/opt/software/multiQC/1.0--singularity/bin/multiqc.img -z -o multiqc/trimmed qc/trimmed"

rule download_ref:
    output:
        "refGenome/ARS-UCD1.2_chr.fa"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        '''
        mkdir refGenome && cd refGenome
        wget --timestamping 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/263/795/GCA_002263795.2_ARS-UCD1.2/GCA_002263795.2_ARS-UCD1.2_genomic.fna.gz' -O ARS-UCD1.2.fa.gz
        gunzip ARS-UCD1.2.fa.gz
        cat ARS-UCD1.2.fa | awk '{if($1 ~ ">CM"){print $0;f=1;}else if($1 ~ ">NKLS"){f=0;}else if(f){print $0;}}' > ARS-UCD1.2_chr.fa
        '''

rule bwa_index:
    input:
        "refGenome/ARS-UCD1.2_chr.fa"
    output:
        expand("refGenome/BwaIndex/genome.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    shell:
        '''
        if [ ! -f refGenome/BwaIndex/genome.fa ];then ln -s ../ARS-UCD1.2_chr.fa refGenome/BwaIndex/genome.fa;fi
        module load bwa/0.7.7.r441
        bwa index -p refGenome/BwaIndex/genome -a bwtsw {input}
        '''

rule bwa_map:
    input:
        expand("refGenome/BwaIndex/genome.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
        r1="data/trimmed/{fragment}_R1_001.fastq.gz",
        r2="data/trimmed/{fragment}_R2_001.fastq.gz",
    output:
#        "{fragment}.test.txt"
        "data/mapped_reads/{fragment}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    threads:4
    shell:
        '''
        name=$(basename {input.r1})
        SM=$(echo $name | cut -d "_" -f1)
        repLib="1"
        if [ -f "data/replicates_list.tab" ];then
            repLib=$(grep $SM $replicates_list | awk '{{ print $2 }}');fi
        LB=$SM.$repLib
        PL="Illumina"
        ##read Fastq 1st read, check the format.If typical, identify ID as "<instrument>:<run number>:<flowcell ID>:<lane>", and PU as the "<instrument>"
        header=$(head -n1 <(zcat {input.r1}) | grep ':*:*:*:*:*:*')
        if [ "$header" != "" ]; then
            PU=$(echo "$header" | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)            ##platform unit-lane ID
        else # "make unique ID and PU using checksum"
            checksum=$(shasum $sample | awk '{{ print $1 }}')
            PU="UnChrPU_"$checksum
        fi
        RGID=$PU.$SM
        echo RGID $RGID LB $LB PL $PL PU $PU SM $SM >> test.txt

        module load bwa/0.7.7.r441
        bwa mem -t {threads} -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" refGenome/BwaIndex/genome {input.r1} {input.r2} > {output}
        '''
