#from os import path
import glob

fastq_pattern = "data/fastq/*/*_001.fastq.gz"
fastq_files = glob.glob(fastq_pattern)
read_ids = [fq[11:-9] for fq in fastq_files]
fragment_ids = list(set([fq[11:-16] for fq in fastq_files]))

Logan = 'AVEAY0010B'.split()
Buri = 'AVEAY001A'.split()
Buri_calf  = 'AVEAY001B AVEAY002A AVEAY002B AVEAY003A AVEAY003B AVEAY004A'.split()
Buri_calf_dam = 'AVEAY004B AVEAY005A AVEAY005B AVEAY006A AVEAY006B AVEAY007A'.split()
HH_calf = 'AVEAY007B AVEAY008A AVEAY008B'.split()
HH_calf_dam = 'AVEAY009A AVEAY009B AVEAY0010A'.split()
HH_calf_sire = 'AVEAY011A AVEAY011B'.split()
Logan_calf = 'AVEAY012A AVEAY012B AVEAY14A'.split()
Logan_calf_dam = 'AVEAY013A AVEAY013B AVEAY14B'.split()

newBuri = 'AVEAY15B'.split()
newBuri_calf_dam = 'AVEAY15A AVEAY16A'.split()

#SAMPLES = Logan + Buri + Buri_calf + Buri_calf_dam + HH_calf + HH_calf_dam
#ID = ['1', '2', '3', '4', '5']
#old_samples = 0
#scale = 8

#SAMPLES = HH_calf_sire + Logan_calf + Logan_calf_dam
#ID = ['6', '7']
#old_samples = 5
#scale = 8

SAMPLES = Logan + Buri + Buri_calf + Buri_calf_dam + HH_calf + HH_calf_dam + HH_calf_sire + Logan_calf + Logan_calf_dam + newBuri + newBuri_calf_dam
ID = ['1', '2', '3', '4', '5', '6', '7', '8']
old_samples = 0
scale = 8

ID2 = ['1', '2']
old_samples2 = 0
scale2 = 8

rule all:
    input:
        expand("qc/fastqc/fastq/{read}_fastqc.{ext}", read=read_ids, ext=['html', 'zip']),
        "qc/multiqc/fastq/multiqc_report.html",
        expand('data/trimmed/{fragment}_{R}_001.fastq.gz', fragment=fragment_ids, R=['R1', 'R2']),
        expand("qc/fastqc/trimmed/{read}_fastqc.{ext}", read=read_ids, ext=['html', 'zip']),
        "qc/multiqc/trimmed/multiqc_report.html",
        "refGenome/ARS-UCD1.2_chr.fa",
        expand("refGenome/BwaIndex/genome.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
        #expand("{fragment}.test.txt", fragment=fragment_ids),
        expand("data/mapped_reads/{fragment}.bam", fragment=fragment_ids),
        expand("data/sorted_reads/{fragment}.bam", fragment=fragment_ids),
        expand("data/merged_reps/{sample}.bam", sample=SAMPLES),
        expand("qc/mappingQC/{sample}.cov", sample=SAMPLES), expand("qc/mappingQC/{sample}.stat", sample=SAMPLES),
        expand("data/dedup/{sample}.bam", sample=SAMPLES), expand("data/dedup/{sample}.txt", sample=SAMPLES),
        expand("data/dedup/{sample}.bai", sample=SAMPLES),
        "refGenome/gatkIndex/genome.fa.fai", "refGenome/gatkIndex/genome.dict",
        "knowVar/bos_taurus_SNPs.vcf", "knowVar/bos_taurus_indels.vcf",
        expand("data/recalib/{sample}.txt", sample=SAMPLES),
        expand("data/recalib/{sample}.bam", sample=SAMPLES),
        expand("vc/hapCaller_single/{sample}.g.vcf", sample=SAMPLES),
        expand("vc/combineGVCF/combined_patch_{id}.vcf", id=ID),
        expand("vc/combineGVCF2/combined_patch_{id2}.vcf", id2=ID2),
        "vc/combined_patch.g.vcf",
        "vc/hapCaller_raw.vcf",
        "vc/hapCaller_raw.SNPs.vcf",
        "vc/hapCaller_raw.INDELs.vcf",
        expand("vc/hapCaller_filtered.{var}s.vcf", var=['SNP','INDEL']),
        expand("vc/hapCaller_pass.{var}s.vcf", var=['SNP','INDEL']),
        "vc/hapCaller_pass_INDELs.monoAllel_edit.vcf",
##        "vc/hapCaller_raw_withExternal.vcf",
##        "vc/hapCaller_raw_withExternal.raw_SNPs.vcf",
##        "vc/hapCaller_raw_withExternal.raw_INDELs.vcf",
##        "vc/hapCaller_raw_withExternal.filtered_SNPs.vcf",
##        "vc/hapCaller_raw_withExternal.filtered_INDELs.vcf",
##        expand("vc/hapCaller_raw_withExternal.filtered_{var}s.vcf", var=['SNP','INDEL']),
##        expand("vc/hapCaller_raw_withExternal.pass_{var}s.vcf", var=['SNP','INDEL']),
##        "vc/hapCaller_raw_withExternal.pass_INDELs.monoAllel_edit.vcf",
        expand("targetReads_grep/singleFiles/{read}_target_reads.fa", read=read_ids),

rule fastqc_pre:
    input:
        "data/fastq/{read}.fastq.gz"
    output:
        html="qc/fastqc/fastq/{read}_fastqc.html",
        zip="qc/fastqc/fastq/{read}_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    threads:1
    shell:
        '''
        #module load FastQC/0.11.5
        module load FastQC/0.11.7-Java-1.8.0_162
        fastqc -t {threads} -f fastq -noextract -o qc/fastqc/fastq/$(dirname {wildcards.read}) {input}
        '''

rule multiQC_pre:
    input:
        html=expand("qc/fastqc/fastq/{read}_fastqc.html", read=read_ids)
    output:
        "qc/multiqc/fastq/multiqc_report.html",
        "qc/multiqc/fastq/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    shell:
        '''
        source activate multiQC
        multiqc -z -o qc/multiqc/fastq qc/fastqc/fastq
        '''
        #"/opt/software/multiQC/1.0--singularity/bin/multiqc.img -z -o qc/multiqc/fastq qc/fastqc/fastq"

rule trimmomatic_pe:
    input:
        r1="data/fastq/{fragment}_R1_001.fastq.gz",
        r2="data/fastq/{fragment}_R2_001.fastq.gz",
    output:
        r1="data/trimmed/{fragment}_R1_001.fastq.gz",
        r2="data/trimmed/{fragment}_R2_001.fastq.gz",
        sr1="data/trimmed/{fragment}_R1.unpaired_001.fastq.gz",
        sr2="data/trimmed/{fragment}_R2.unpaired_001.fastq.gz",
#    log:
#        "logs/trimmomatic/{fragment}.log"
#    params:
#        # list of trimmers (see manual)
#        trimmer=["TRAILING:3"],
#        # optional parameters
#        extra=""
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 16
    threads:2
    shell:
        '''
        #module load Trimmomatic/0.33
        #java -jar $TRIM/trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} {output.r1} {output.sr1} {output.r2} {output.sr2} ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:20
        module load Trimmomatic/0.38-Java-1.8.0_162
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads {threads} -phred33 {input.r1} {input.r2} {output.r1} {output.sr1} {output.r2} {output.sr2} ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:20
        '''

rule fastqc_post:
    input:
        "data/trimmed/{read}.fastq.gz"
    output:
        html="qc/fastqc/trimmed/{read}_fastqc.html",
        zip="qc/fastqc/trimmed/{read}_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    threads:1
    shell:
        '''
        #module load FastQC/0.11.5
        module load FastQC/0.11.7-Java-1.8.0_162
        fastqc -t {threads} -f fastq -noextract -o qc/fastqc/trimmed/$(dirname {wildcards.read}) {input}
        '''

rule multiQC_post:
    input:
        html=expand("qc/fastqc/trimmed/{read}_fastqc.html", read=read_ids)
    output:
        "qc/multiqc/trimmed/multiqc_report.html",
        "qc/multiqc/trimmed/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    shell:
        '''
        source activate multiQC
        multiqc -z -o qc/multiqc/trimmed qc/fastqc/trimmed
        '''
        #"/opt/software/multiQC/1.0--singularity/bin/multiqc.img -z -o qc/multiqc/trimmed qc/fastqc/trimmed"

rule download_ref:
    output:
        "refGenome/ARS-UCD1.2_chr.fa"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
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
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 16
    shell:
        '''
        if [ ! -f refGenome/BwaIndex/genome.fa ];then ln -s ../ARS-UCD1.2_chr.fa refGenome/BwaIndex/genome.fa;fi
        #module load bwa/0.7.7.r441
        module swap OpenMPI  OpenMPI/2.1.1
        module load BWA/0.7.17
        bwa index -p refGenome/BwaIndex/genome -a bwtsw {input}
        '''

## This rule fail to add the appropriate RG info (only for the samples in 10_*). I discovered this after samtools_sort step. check "correct_RG.sh" 
rule bwa_map:
    input:
        expand("refGenome/BwaIndex/genome.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
        r1="data/trimmed/{fragment}_R1_001.fastq.gz",
        r2="data/trimmed/{fragment}_R2_001.fastq.gz",
    output:
        #"{fragment}.test.txt"
        "data/mapped_reads/{fragment}.bam",   ## Actually the output is sam file
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 24,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 32
    threads:4
    shell:
        '''
        name=$(basename {input.r1})
        SM=$(echo $name | cut -d "_" -f1)
        LB=$(echo $name | cut -d"_" -f1,2)  ## We use <Index.Sequence> in the Illumina file name as an index to the library
        batch=$(basename "$(dirname {input.r1})")
        if [ "$batch" != "trimmed" ];then LB=$batch.$LB;fi
        #repLib="1"
        #if [ -f "data/replicates_list.tab" ];then
        #    repLib=$(grep $SM $replicates_list | awk '{{ print $2 }}');fi
        #LB=$SM.$repLib  ## This way we use metadeta file to track library info. e.g. you have multiple libraries with the same index
        PL="Illumina"
        ##read Fastq 1st read, check the format.If typical, identify ID as "<instrument>:<run number>:<flowcell ID>:<lane>"
        header=$(head -n1 <(zcat {input.r1}) | grep ':*:*:*:*:*:*')
        #if [ "$header" != "" ]; then
        #    PU=$(echo "$header" | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)       
        #else # "make unique ID and PU using checksum"
        #    checksum=$(shasum $sample | awk '{{ print $1 }}')
        #    PU="UnChrPU_"$checksum
        #fi
        #RGID=$PU.$SM
        if [ "$header" != "" ]; then
            RGID=$(echo "$header" | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)     
        else # "make unique ID and PU using checksum"
            checksum=$(shasum {input.r1} | awk '{{ print $1 }}')
            RGID="UnChrPU_"$checksum
        fi
        PU=$RGID.$LB  
        #echo RGID $RGID LB $LB PL $PL PU $PU SM $SM > {output}
        echo RGID $RGID LB $LB PL $PL PU $PU SM $SM >> test.txt

        #module load bwa/0.7.7.r441
        module swap OpenMPI  OpenMPI/2.1.1
        module load BWA/0.7.17
        #bwa mem -t {threads} -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" refGenome/BwaIndex/genome {input.r1} {input.r2} > {output}
        bwa mem -t {threads} -M -R '@RG\\tID:$RGID\\tSM:$SM\\tPL:$PL\\tLB:$LB\\tPU:$PU' refGenome/BwaIndex/genome {input.r1} {input.r2} > {output} 
        '''

rule samtools_sort:
    input:
        "data/mapped_reads/{fragment}.bam",
    output:
        "data/sorted_reads/{fragment}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 8,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 64
    threads:8
    shell:
        '''
        rm -f data/sorted_reads/{wildcards.fragment}*
        module swap GNU GNU/5.4.0-2.26
        module swap OpenMPI OpenMPI/1.10.3
        module load SAMtools/1.5
        #module load SAMTools/1.5
        samtools sort -T data/sorted_reads/{wildcards.fragment} -m 4G -@ {threads} -O bam {input} > {output}
        '''

# There is a problem with bwa which add the mapping command to the @PG line. Their current implementation translate \t into real tabs ending up with 2 ID keys in one @PG line
# https://github.com/samtools/htsjdk/issues/677 
# https://github.com/samtools/hts-specs/issues/275
# https://github.com/broadinstitute/picard/issues/1139
# To work around this problem, I use "VALIDATION_STRINGENCY=LENIENT" in picardTools commands
#
## shortcut code to avoid schduling jobs for single files
# cd Bovine_seq
### for sample in AVEAY011A AVEAY011B AVEAY012A AVEAY012B AVEAY013A AVEAY013B AVEAY14A AVEAY14B;do
# for sample in AVEAY14A AVEAY14B;do
#   ln -s $(pwd)/data/sorted_reads/*/${sample}_*.bam data/merged_reps/${sample}.bam
# done
#
rule merge_rep:
#    input:
#        expand("data/sorted_reads/{fragment}.bam", fragment=fragment_ids),
    output:
        "data/merged_reps/{sample}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    shell:
        '''
        inputs_bam=()
        while read bam;do
            input_bam=$"INPUT="$bam;
            inputs_bam+=($input_bam);
        done < <(find data/sorted_reads -name "{wildcards.sample}_*")
        if [ ${{#inputs_bam[@]}} -eq 1 ];then
            singleSample=${{inputs_bam[0]#INPUT=}}
            ln -s $(pwd)/$singleSample {output}
        else
            module load picard/2.18.1-Java-1.8.0_152
            java -Xmx8g -jar $EBROOTPICARD/picard.jar MergeSamFiles VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true ${{inputs_bam[*]}} OUTPUT={output}
        fi
        '''

rule mappingQC:
    input:
        "data/merged_reps/{sample}.bam",
    output:
        cov="qc/mappingQC/{sample}.cov",
        stat="qc/mappingQC/{sample}.stat",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    shell:
        '''
        module swap GNU GNU/5.4.0-2.26
        module swap OpenMPI OpenMPI/1.10.3
        module load SAMtools/1.5
        samtools depth {input} | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' > {output.cov}
        samtools flagstat {input} > {output.stat}
        '''

rule markDuplicates:
    input:
        "data/merged_reps/{sample}.bam",
    output:
        bam="data/dedup/{sample}.bam",
        metrics="data/dedup/{sample}.txt",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 8,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 12
    shell:
        '''
        module load picard/2.18.1-Java-1.8.0_152
        java -Xmx8g -jar $EBROOTPICARD/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}
        '''

rule buildBamIndex:
    input:
        "data/dedup/{sample}.bam",
    output:
        "data/dedup/{sample}.bai",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    shell:
        '''
        module load picard/2.18.1-Java-1.8.0_152
        java -Xmx8g -jar $EBROOTPICARD/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT={input}
        '''

rule GATK_index:
    input:
        "refGenome/ARS-UCD1.2_chr.fa"
    output:
        index="refGenome/gatkIndex/genome.fa.fai",
        dict="refGenome/gatkIndex/genome.dict",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    shell:
        '''
        mkdir -p refGenome/gatkIndex
        if [ ! -f refGenome/gatkIndex/genome.fa ];then ln -s ../ARS-UCD1.2_chr.fa refGenome/gatkIndex/genome.fa;fi
        module swap GNU GNU/5.4.0-2.26
        module swap OpenMPI OpenMPI/1.10.3
        module load SAMtools/1.5
        module load picard/2.18.1-Java-1.8.0_152
        #module load picardTools/1.89
        samtools faidx "refGenome/gatkIndex/genome.fa"
        java -Xmx4g -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R= {input} O= {output.dict}
        #java -Xmx4g -jar $PICARD/CreateSequenceDictionary.jar R= {input} O= {output.dict}
        '''

rule download_knowVar:
    output:
        knownSNPs="knowVar/bos_taurus_SNPs.vcf",
        knownIndels="knowVar/bos_taurus_indels.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    shell:
        '''
        mkdir knowVar && cd knowVar
        wget 'ftp://ftp.ensembl.org/pub/release-94/variation/vcf/bos_taurus/bos_taurus.vcf.gz' -O bos_taurus.vcf.gz
        gunzip bos_taurus.vcf.gz
        cd ../
        grep "^>" refGenome/gatkIndex/genome.fa | awk -F"[>, ]" '{if($15 != "")print $15,$2;else print "MT",$2}' > knowVar/chr_map
        module load GCC/4.9.3-2.25
        module swap OpenMPI OpenMPI/1.10.2
        module load tabix/0.2.6
        cat knowVar/bos_taurus.vcf | bgzip > knowVar/bos_taurus.vcf.gz
        tabix -p vcf knowVar/bos_taurus.vcf.gz
        module swap GCC GCC/6.4.0-2.28
        module swap OpenMPI OpenMPI/2.1.2
        module load bcftools/1.9.64
        bcftools annotate --rename-chrs knowVar/chr_map -o knowVar/bos_taurus_fixedChrNames.vcf knowVar/bos_taurus.vcf.gz
        perl scripts/sortByRef.pl knowVar/bos_taurus_fixedChrNames.vcf refGenome/gatkIndex/genome.fa.fai > knowVar/bos_taurus_fixedChrNames_sorted.vcf
        grep "^#" knowVar/bos_taurus.vcf > knowVar/bos_taurus_SNPs.vcf
        grep "TSA=SNV" knowVar/bos_taurus_fixedChrNames_sorted.vcf >> knowVar/bos_taurus_SNPs.vcf
        grep "^#" knowVar/bos_taurus.vcf > knowVar/bos_taurus_indels.vcf
        grep -v "TSA=SNV" knowVar/bos_taurus_fixedChrNames_sorted.vcf >> knowVar/bos_taurus_indels.vcf
        module load Java/1.8.0_172
        source activate gatk
        gatk IndexFeatureFile -F knowVar/bos_taurus_SNPs.vcf
        gatk IndexFeatureFile -F knowVar/bos_taurus_indels.vcf
        '''

rule BaseRecalib:
    input:
        bam="data/dedup/{sample}.bam",
        ref="refGenome/gatkIndex/genome.fa",
        knownSNPs="knowVar/bos_taurus_SNPs.vcf",
        knownIndels="knowVar/bos_taurus_indels.vcf",
    output:
        report="data/recalib/{sample}.txt",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 16
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx15G" BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        --use-original-qualities \
        -O {output.report} \
        --known-sites {input.knownSNPs} \
        --known-sites {input.knownIndels}
        '''

rule ApplyBQSR:
    input:
        bam="data/dedup/{sample}.bam",
        ref="refGenome/gatkIndex/genome.fa",
        report="data/recalib/{sample}.txt",
    output:
        bam="data/recalib/{sample}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 10,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 16
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx15G" ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.bam} \
        -bqsr {input.report} \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
        --add-output-sam-program-record \
        --use-original-qualities
        '''

# https://software.broadinstitute.org/gatk/documentation/article?id=3893
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.9.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php
# https://software.broadinstitute.org/gatk/documentation/article?id=11073
# https://software.broadinstitute.org/gatk/documentation/article?id=11068
rule HaplotypeCaller_single:
    input:
        "refGenome/gatkIndex/genome.fa.fai", "refGenome/gatkIndex/genome.dict",
        ref="refGenome/gatkIndex/genome.fa",
        bam="data/recalib/{sample}.bam",
    output:
        gvcf="vc/hapCaller_single/{sample}.g.vcf",
        idx="vc/hapCaller_single/{sample}.g.vcf.idx",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 168,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 32
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx30G" HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        --emit-ref-confidence GVCF \
        --pcr-indel-model NONE \
        -O {output.gvcf}
        '''

## if this is not the first batch, select the appropriate SAMPLES
rule combineGVCFs:
    input:
        gvcfs=expand("vc/hapCaller_single/{sample}.g.vcf", sample=SAMPLES),
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/combineGVCF/combined_patch_{id}.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 24 * 3,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 64
    threads:1
    shell:
        '''
        inputs=()
        for gvcf in {input.gvcfs};do
            input_gvcf=$" -V "$gvcf;
            inputs+=($input_gvcf);
        done

        start_id=$((({wildcards.id}-{old_samples}-1)*{scale}))
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx62G" CombineGVCFs \
        -R {input.ref} \
        ${{inputs[*]:$start_id:{scale}}} \
        -O {output}
        '''

rule combineGVCFs2:
    input:
        gvcfs=expand("vc/combineGVCF/combined_patch_{id}.vcf", id=ID),
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/combineGVCF2/combined_patch_{id2}.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 24 * 3,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 64
    threads:1
    shell:
        '''
        inputs=()
        for gvcf in {input.gvcfs};do
            input_gvcf=$" -V "$gvcf;
            inputs+=($input_gvcf);
        done

        start_id=$((({wildcards.id2}-{old_samples2}-1)*{scale2}))
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx62G" CombineGVCFs \
        -R {input.ref} \
        ${{inputs[*]:$start_id:{scale2}}} \
        -O {output}
        '''

rule combineGVCFs_final:
    input:
        gvcfs=expand("vc/combineGVCF2/combined_patch_{id2}.vcf", id2=ID2),
        ext_gvcfs=expand("external/SRP072240/vc/combineGVCF/combined_patch_{id}.vcf", id=['1']),
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/combined_patch.g.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 24 * 5,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 64
    threads:1
    shell:
        '''
        inputs=()
        for gvcf in {input.gvcfs} {input.ext_gvcfs};do
            input_gvcf=$" -V "$gvcf;
            inputs+=($input_gvcf);
        done
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx62G" CombineGVCFs \
        -R {input.ref} \
        ${{inputs[*]}} \
        -O {output}
        '''

## if this is not the first batch, revert to the whole list SAMPLES and run the touch script to avoid re-runs
rule genotypeGVCFs:
    input:
        gvcf="vc/combined_patch.g.vcf",
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/hapCaller_raw.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 24 * 7,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 64
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx60G" GenotypeGVCFs \
        -R {input.ref} \
        -V {input.gvcf} \
        --max-alternate-alleles 6 \
        -O {output}
        '''

rule select_snps:
    input:
        ref="refGenome/gatkIndex/genome.fa",
        vcf="vc/hapCaller_raw.vcf",
    output:
        "vc/hapCaller_raw.SNPs.vcf"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 24
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx60G" SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type "SNP" \
        -O {output}
        '''

rule select_indels:
    input:
        ref="refGenome/gatkIndex/genome.fa",
        vcf="vc/hapCaller_raw.vcf",
    output:
        "vc/hapCaller_raw.INDELs.vcf"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 24
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx60G" SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type "INDEL" \
        -O {output}
        '''

rule filter_variants:
    input:
        expand("vc/hapCaller_raw.{var}s.vcf", var=['SNP','INDEL']),
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/hapCaller_filtered.{var}s.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 24
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx20G" VariantFiltration \
        -R {input.ref} \
        -V "vc/hapCaller_raw.{wildcards.var}s.vcf" \
        --filter-name "filterQD" --filter-expression "vc.hasAttribute('QD') && QD < 3.0" \
        --filter-name "filterFS" --filter-expression "vc.hasAttribute('FS') && FS > 50.0" \
        --filter-name "filterSOR" --filter-expression "vc.hasAttribute('SOR') && SOR > 3.0" \
        --filter-name "filterMQ" --filter-expression "vc.hasAttribute('MQ') && MQ < 50.0" \
        --filter-name "filterMQRankSum" --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -2.0" \
        --filter-name "filterHiReadPosRankSum" --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -4.0" \
        --filter-name "filterLowReadPosRankSum" --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum > 5.0" \
        -O {output}
        '''

rule apply_filters:
    input:
        expand("vc/hapCaller_filtered.{var}s.vcf", var=['SNP','INDEL']),
    output:
        "vc/hapCaller_pass.{var}s.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 24
    threads:1
    shell:
        '''
        grep -E "^#|PASS" "vc/hapCaller_filtered.{wildcards.var}s.vcf" > "vc/hapCaller_pass.{wildcards.var}s.vcf"
        grep -v "PASS" "vc/hapCaller_filtered.{wildcards.var}s.vcf" > "vc/hapCaller_failed.{wildcards.var}s.vcf"
        '''


rule prep_indels:
    input:
        "vc/hapCaller_pass.INDELs.vcf"
    output:
        monoAllel_vcf="vc/hapCaller_pass.INDELs.monoAllel.vcf",
        edit_vcf="vc/hapCaller_pass_INDELs.monoAllel_edit.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 24
    threads:1
    shell:
        '''
        awk '/#/{{print;next}}{{if($5 !~ /,/){{print}}}}' {input} > {output.monoAllel_vcf}
        awk 'BEGIN{{FS="\t";OFS="\t"}}/#/{{print;next}}{{if(length($4)>1){{$5!="A"?$4="A":$4="T";}};if(length($5)>1){{$4!="A"?$5="A":$5="T";}};print;}}' {output.monoAllel_vcf} > {output.edit_vcf}
        '''

rule findTargetReads:
    input:
        seq="data/trimmed/{read}.fastq.gz",
        bait="targetReads_grep/bait_Allkmers.txt",
    output:
        fq="targetReads_grep/singleFiles/{read}_target_reads.fq",
        fa="targetReads_grep/singleFiles/{read}_target_reads.fa",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    threads:1
    shell:
        '''
        zgrep -B1 -A2 -F -f {input.bait} {input.seq} | grep -v "^\-\-" > {output.fq} || true
        cat {output.fq} | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > {output.fa}
        '''

