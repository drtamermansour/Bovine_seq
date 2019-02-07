#from os import path
import glob

Logan = 'AVEAY0010B'.split()
Buri = 'AVEAY001A'.split()
Buri_calf  = 'AVEAY001B AVEAY002A AVEAY002B AVEAY003A AVEAY003B AVEAY004A'.split()
Buri_calf_dam = 'AVEAY004B AVEAY005A AVEAY005B AVEAY006A AVEAY006B AVEAY007A'.split()
HH_calf = 'AVEAY007B AVEAY008A AVEAY008B'.split()
HH_calf_dam = 'AVEAY009A AVEAY009B AVEAY0010A'.split()
HH_calf_sire = 'AVEAY011A AVEAY011B'.split()
Logan_calf = 'AVEAY012A AVEAY012B AVEAY14A'.split()
Logan_calf_dam = 'AVEAY013A AVEAY013B AVEAY14B'.split()

#SAMPLES = Logan + Buri + Buri_calf + Buri_calf_dam + HH_calf + HH_calf_dam
#ID = ['1', '2', '3', '4', '5']
#old_samples = 0
#scale = 8

#SAMPLES = HH_calf_sire + Logan_calf + Logan_calf_dam
#ID = ['6', '7']
#old_samples = 5
#scale = 8

SAMPLES = Logan + Buri + Buri_calf + Buri_calf_dam + HH_calf + HH_calf_dam + HH_calf_sire + Logan_calf + Logan_calf_dam
ID = ['1', '2', '3', '4', '5', '6', '7']

fastq_pattern = "data/fastq/*/*_001.fastq.gz"
fastq_files = glob.glob(fastq_pattern)
read_ids = [fq[11:-9] for fq in fastq_files]
fragment_ids = list(set([fq[11:-16] for fq in fastq_files]))

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
        #expand("data/merged_reps/{sample}.bai", sample=SAMPLES),
        expand("data/dedup/{sample}.bai", sample=SAMPLES),
        "refGenome/gatkIndex/genome.fa.fai", "refGenome/gatkIndex/genome.dict",
        expand("vc/hapCaller_single/{sample}.g.vcf", sample=SAMPLES),
        expand("vc/combineGVCF/combined_patch_{id}.vcf", id=ID),
#        "vc/hapCaller_raw.vcf",
        "vc/hapCaller_raw_withExternal.vcf",
        "vc/hapCaller_raw_withExternal.raw_SNPs.vcf",
        "vc/hapCaller_raw_withExternal.raw_INDELs.vcf",
        "vc/hapCaller_raw_withExternal.filtered_SNPs.vcf",
        "vc/hapCaller_raw_withExternal.filtered_INDELs.vcf",
        expand("vc/hapCaller_raw_withExternal.filtered_{var}s.vcf", var=['SNP','INDEL']),
        expand("vc/hapCaller_raw_withExternal.pass_{var}s.vcf", var=['SNP','INDEL']),
        "vc/hapCaller_raw_withExternal.pass_INDELs.monoAllel_edit.vcf",
        expand("targetReads_grep/singleFiles/{read}_target_reads.fa", read=read_ids),

rule fastqc_pre:
    input:
        "data/fastq/{read}.fastq.gz"
    output:
        html="qc/fastqc/fastq/{read}_fastqc.html",
        zip="qc/fastqc/fastq/{read}_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    threads:1
    shell:
        '''
        module load FastQC/0.11.5
        fastqc -t {threads} -f fastq -noextract -o qc/fastqc/fastq/$(dirname {wildcards.read}) {input}
        '''

rule multiQC_pre:
    input:
        html=expand("qc/fastqc/fastq/{read}_fastqc.html", read=read_ids)
    output:
        "qc/multiqc/fastq/multiqc_report.html",
        "qc/multiqc/fastq/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        "/opt/software/multiQC/1.0--singularity/bin/multiqc.img -z -o qc/multiqc/fastq qc/fastqc/fastq" 

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
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 8,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 4
    threads:2
    shell:
        '''
        module load Trimmomatic/0.33
        java -jar $TRIM/trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} {output.r1} {output.sr1} {output.r2} {output.sr2} ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:20
        '''

rule fastqc_post:
    input:
        "data/trimmed/{read}.fastq.gz"
    output:
        html="qc/fastqc/trimmed/{read}_fastqc.html",
        zip="qc/fastqc/trimmed/{read}_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    threads:1
    shell:
        '''
        module load FastQC/0.11.5
        fastqc -t {threads} -f fastq -noextract -o qc/fastqc/trimmed/$(dirname {wildcards.read}) {input}
        '''

rule multiQC_post:
    input:
        html=expand("qc/fastqc/trimmed/{read}_fastqc.html", read=read_ids)
    output:
        "qc/multiqc/trimmed/multiqc_report.html",
        "qc/multiqc/trimmed/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        "/opt/software/multiQC/1.0--singularity/bin/multiqc.img -z -o qc/multiqc/trimmed qc/fastqc/trimmed"

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
        #"{fragment}.test.txt"
        "data/mapped_reads/{fragment}.bam",   ## Actually the output is sam file
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 24,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 120
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

        module load bwa/0.7.7.r441
        bwa mem -t {threads} -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" refGenome/BwaIndex/genome {input.r1} {input.r2} > {output}
        '''

rule samtools_sort:
    input:
        "data/mapped_reads/{fragment}.bam",
    output:
        "data/sorted_reads/{fragment}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 8,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 64
    threads:4
    shell:
        '''
        rm -f data/sorted_reads/{wildcards.fragment}*
        module load SAMTools/1.5
        samtools sort -T data/sorted_reads/{wildcards.fragment} -m 2G -@ {threads} -O bam {input} > {output}
        '''

# There is a problem with bwa which add the mapping command to the @PG line. Their current implementation translate \t into real tabs ending up with 2 ID keys in one @PG line
# https://github.com/samtools/htsjdk/issues/677 
# https://github.com/samtools/hts-specs/issues/275
# https://github.com/broadinstitute/picard/issues/1139
# To work around this problem, I use "VALIDATION_STRINGENCY=LENIENT" in picardTools commands
#
## shortcut code to avoid schduling jobs for single files
# cd Bovine_seq
# for sample in AVEAY011A AVEAY011B AVEAY012A AVEAY012B AVEAY013A AVEAY013B AVEAY14A AVEAY14B;do
#   ln -s $(pwd)/data/sorted_reads/*/${sample}_*.bam data/merged_reps/${sample}.bam
# done
#
rule merge_rep:
    input:
        expand("data/sorted_reads/{fragment}.bam", fragment=fragment_ids),
    output:
        "data/merged_reps/{sample}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
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
            module load picardTools/1.89
            java -Xmx8g -jar $PICARD/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true ${{inputs_bam[*]}} OUTPUT={output}
        fi
        '''
#    shell:
#        '''
#        module load picardTools/1.89
#        inputs_bam=()
#        while read bam;do
#            input_bam=$"INPUT="$bam;
#            inputs_bam+=($input_bam);
#        done < <(find data/sorted_reads -name "{wildcards.sample}_*")
#        java -Xmx8g -jar $PICARD/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true ${{inputs_bam[*]}} OUTPUT={output}
#        '''

rule mappingQC:
    input:
        "data/merged_reps/{sample}.bam",
    output:
        cov="qc/mappingQC/{sample}.cov",
        stat="qc/mappingQC/{sample}.stat",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    shell:
        '''
        module load SAMTools/1.5
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
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 8,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 12
    shell:
        '''
        module load picardTools/1.89
        java -Xmx8g -jar $PICARD/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}
        '''

rule buildBamIndex:
    input:
        "data/dedup/{sample}.bam",
#        "data/merged_reps/{sample}.bam",
    output:
        "data/dedup/{sample}.bai",
#        "data/merged_reps/{sample}.bai",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    shell:
        '''
        module load picardTools/1.89
        java -Xmx8g -jar $PICARD/BuildBamIndex.jar VALIDATION_STRINGENCY=LENIENT INPUT={input}
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
        if [ ! -f refGenome/gatkIndex/genome.fa ];then ln -s ../ARS-UCD1.2_chr.fa refGenome/gatkIndex/genome.fa;fi
        module load SAMTools/1.5
        module load picardTools/1.89
        samtools faidx "refGenome/gatkIndex/genome.fa"
        java -Xmx4g -jar $PICARD/CreateSequenceDictionary.jar R= {input} O= {output.dict}
        '''

rule HaplotypeCaller_single:
    input:
        "refGenome/gatkIndex/genome.fa.fai", "refGenome/gatkIndex/genome.dict",
        ref="refGenome/gatkIndex/genome.fa",
        bam="data/dedup/{sample}.bam",
    output:
        gvcf="vc/hapCaller_single/{sample}.g.vcf", 
        idx="vc/hapCaller_single/{sample}.g.vcf.idx",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 96,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 64
    threads:4
    shell:
        '''
        module load GATK/3.7.0
        java -Xmx28g -jar $GATK/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        --emitRefConfidence GVCF \
        -nct 3 \
        -o {output.gvcf}
#        --dbsnp $snps \
#        --variant_index_type LINEAR \
#        --variant_index_parameter 128000 \
#        -nct 3 \
#        -o {output}
        '''

## if this is not the first batch, select the appropriate SAMPLES
rule combineGVCFs:
    input:
        gvcfs=expand("vc/hapCaller_single/{sample}.g.vcf", sample=SAMPLES),
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/combineGVCF/combined_patch_{id}.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 24 * 3,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 128
    threads:1
    shell:
        '''
        touch vc/hapCaller_single/*.g.vcf.idx
        inputs=()
        #for gvcf in vc/hapCaller_single/*.g.vcf;do
        for gvcf in {input.gvcfs};do
            input_gvcf=$" -V "$gvcf;
            inputs+=($input_gvcf);
        done

        #start_id=$((({wildcards.id}-1)*8))
        start_id=$((({wildcards.id}-{old_samples}-1)*{scale}))
        module load GATK/3.7.0
        java -Xmx120g -jar $GATK/GenomeAnalysisTK.jar \
        -T CombineGVCFs \
        -R {input.ref} \
        ${{inputs[*]:$start_id:{scale}}} \
        -o {output}
#        $(echo $trim_inputs) \
#        -o {output}
        '''

## if this is not the first batch, revert to the whole list SAMPLES and run the touch script to avoid re-runs
rule genotypeGVCFs:
    input:
        #expand("vc/hapCaller_single/{sample}.g.vcf", sample=SAMPLES),
        expand("vc/combineGVCF/combined_patch_{id}.vcf", id=ID),
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/hapCaller_raw.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 24 * 5,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 128
    threads:1
    shell:
        '''
        #touch vc/hapCaller_single/*.g.vcf.idx
        touch vc/combineGVCF/*.vcf.idx
        inputs=()
        #for gvcf in vc/hapCaller_single/*.g.vcf;do
        for gvcf in vc/combineGVCF/combined_patch_*.vcf;do
            input_gvcf=$" -V "$gvcf;
            inputs+=($input_gvcf);
        done

        module load GATK/3.7.0
        java -Xmx120g -jar $GATK/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R {input.ref} \
        ${{inputs[*]}} \
        --max_alternate_alleles 6 \
        -o {output}
#        -nt 3 \     								## multithreading cause runtime error
#        --disable_auto_index_creation_and_locking_when_reading_rods \          ## no need because I am touching the indices 
#        -o {output}
#        --dbsnp $snps \
#        $(echo $samples) \
        '''

## run the Snakmeke files of external experiments before you proceed
rule genotypeGVCFs2:
    input:
        expand("vc/combineGVCF/combined_patch_{id}.vcf", id=ID),
        expand("external/SRP072240/vc/combineGVCF/combined_patch_{id}.vcf", id=['1']),
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/hapCaller_raw_withExternal.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 24 * 5,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 128
    threads:1
    shell:
        '''
        #touch vc/combineGVCF/*.vcf.idx
        touch external/SRP072240/vc/combineGVCF/*.vcf.idx
        inputs=()
        for gvcf in vc/combineGVCF/combined_patch_*.vcf external/SRP072240/vc/combineGVCF/combined_patch_*.vcf;do
            input_gvcf=$" -V "$gvcf;
            inputs+=($input_gvcf);
        done

        module load GATK/3.7.0
        java -Xmx120g -jar $GATK/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R {input.ref} \
        ${{inputs[*]}} \
        --max_alternate_alleles 6 \
        -o {output}
        '''

rule select_snps:
    input:
        ref="refGenome/gatkIndex/genome.fa",
        vcf="vc/hapCaller_raw_withExternal.vcf",
    output:
        "vc/hapCaller_raw_withExternal.raw_SNPs.vcf"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 24
    threads:4
    shell:
        '''
        module load GATK/3.7.0
        java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -selectType "SNP" \
        -nt 3 \
        -o {output}
        '''

rule select_indels:
    input:
        ref="refGenome/gatkIndex/genome.fa",
        vcf="vc/hapCaller_raw_withExternal.vcf",
    output:
        "vc/hapCaller_raw_withExternal.raw_INDELs.vcf"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 24
    threads:4
    shell:
        '''
        module load GATK/3.7.0
        java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -selectType "INDEL" \
        -nt 3 \
        -o {output}
        '''

rule filter_variants:
    input:
        expand("vc/hapCaller_raw_withExternal.raw_{var}s.vcf", var=['SNP','INDEL']),
        ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/hapCaller_raw_withExternal.filtered_{var}s.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 24
    threads:1
    shell:
        '''
        module load GATK/3.7.0
        java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R {input.ref} \
        -V "vc/hapCaller_raw_withExternal.raw_{wildcards.var}s.vcf" \
        --filterName "filterQD" --filterExpression "vc.hasAttribute('QD') && QD < 3.0" \
        --filterName "filterFS" --filterExpression "vc.hasAttribute('FS') && FS > 50.0" \
        --filterName "filterSOR" --filterExpression "vc.hasAttribute('SOR') && SOR > 3.0" \
        --filterName "filterMQ" --filterExpression "vc.hasAttribute('MQ') && MQ < 50.0" \
        --filterName "filterMQRankSum" --filterExpression "vc.hasAttribute('MQRankSum') && MQRankSum < -2.0" \
        --filterName "filterHiReadPosRankSum" --filterExpression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -4.0" \
        --filterName "filterLowReadPosRankSum" --filterExpression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum > 5.0" \
        -o {output}
        '''

rule apply_filters:
    input:
        expand("vc/hapCaller_raw_withExternal.filtered_{var}s.vcf", var=['SNP','INDEL']),
    output:
        "vc/hapCaller_raw_withExternal.pass_{var}s.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 24
    threads:1
    shell:
        '''
        grep -E "^#|PASS" "vc/hapCaller_raw_withExternal.filtered_{wildcards.var}s.vcf" > "vc/hapCaller_raw_withExternal.pass_{wildcards.var}s.vcf"
        #grep "PASS" "vc/hapCaller_raw_withExternal.raw_{wildcards.var}s.vcf" >> "vc/hapCaller_raw_withExternal.pass_{wildcards.var}s.vcf"
        grep -v "PASS" "vc/hapCaller_raw_withExternal.filtered_{wildcards.var}s.vcf" > "vc/hapCaller_raw_withExternal.failed_{wildcards.var}s.vcf"
        '''

rule prep_indels:
    input:
        "vc/hapCaller_raw_withExternal.pass_INDELs.vcf"
    output:
       monoAllel_vcf="vc/hapCaller_raw_withExternal.pass_INDELs.monoAllel.vcf",
       edit_vcf="vc/hapCaller_raw_withExternal.pass_INDELs.monoAllel_edit.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 24
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
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    threads:1
    shell:
        '''
        zgrep -B1 -A2 -F -f {input.bait} {input.seq} | grep -v "^\-\-" > {output.fq} || true
        cat {output.fq} | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > {output.fa}
        '''

