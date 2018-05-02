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

rule all:
    input:
#        expand('output/{sample}.txt', sample=SAMPLES)#, R=['R1', 'R2'])
        expand("qc/{sample}_fastqc.html", sample=sample_ids),
        "multiqc/multiqc_report.html"

#rule test:
#    input:
#        lambda wildcards: [f for f in fastq_files if wildcards.sample in f]
##        fastq_files
##        "data/fastq/{batch}/{sample}_{index}_{lane}_{R}_001.fastq.gz"
#    output:
##        'output/{sample}_{R}.txt'
#        'output/{sample}.txt'
#    shell:
#        "echo {input} > {output}"


rule fastqc:
    input:
        "data/fastq/{sample}.fastq.gz"
    output:
        dir="qc",
        html="qc/{sample}_fastqc.html",
        zip="qc/{sample}_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    threads:1
    shell:
        '''
        module load FastQC/0.11.5
        fastqc -t {threads} -f fastq -noextract -o qc/$(dirname {wildcards.sample}) {input}
        '''

rule multiQC:
    input:
        dir="qc",
    output:
        "multiqc/multiqc_report.html",
        "multiqc/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        "/opt/software/multiQC/1.0--singularity/bin/multiqc.img -z -o multiqc {input.dir}" 

