work_dir=$(pwd)
mkdir -p ./targetReads_grep_plasmidClone
## data processing
# ....

## find target reads
cd $work_dir/targetReads_grep_plasmidClone ## The folder where you have the bait.fa files

## unwrap the bait sequences 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < bait.fa | tail -n +2 > bait_unwrapped.fa

## convert the bait to kmers
k=25
while read header; read line;do
  step=1;
  end=$((${#line} +1 - $k));
  for i in $(seq 1 $step $end);do
    #echo $step $end $i;
    kmer=`expr substr $line $i $k`;
    echo $kmer;done
    #echo -e "$header.position.$i\n$kmer";done
done < bait_unwrapped.fa | grep -iv n | tr "atcg" "ATCG" | sort | uniq > bait_kmers.txt

## find reverse complement kmers
cat bait_kmers.txt | rev | tr "ATGC" "TACG" > bait_rckmers.txt

## make one list of all possible kmers
cat bait_kmers.txt bait_rckmers.txt | sort | uniq > bait_Allkmers.txt

## find matching reads
cd $work_dir
mkdir -p targetReads_grep_plasmidClone/singleFiles/{10_17_18,10_9_18,12_19_17,1_4_18-1,1_4_18-2,5_2_18,6_13_18,8_2_18}
while read name;do
  read=${name%%.fastq.gz};
  echo $read
  seq="data/trimmed/${read}.fastq.gz"
  bait="targetReads_grep_plasmidClone/bait_Allkmers.txt"
  fq="targetReads_grep_plasmidClone/singleFiles/${read}_target_reads.fq"
  fa="targetReads_grep_plasmidClone/singleFiles/${read}_target_reads.fa"
  #zgrep -B1 -A2 -F -f ${bait} ${seq} | grep -v "^\-\-" > ${fq}
  #cat ${fq} | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > ${fa}
  zgrep -B1 -F -f ${bait} ${seq} | grep -v "^\-\-" | sed 's/^@/>/' > ${fa}
#done < <(ls data/fastq/*/*.fastq.gz | sed 's|data/fastq/||')
done < <(ls data/fastq/*/AVEAY15B_*.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/*/AVEAY003B_*.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/*/AVEAY003A_*.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/10_17_18/AVEAY15B_S14_L006_R2_001.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/*/AVEAY001A_*.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/1_4_18-1/AVEAY001B_S2_L007_R2_001.fastq.gz | sed 's|data/fastq/||')

## recognize the ends
for f in targetReads_grep_plasmidClone/singleFiles/*/*_R1*_*.fa;do sed -i 's|^\(>.*\) |\1\\1 |' $f;done
for f in targetReads_grep_plasmidClone/singleFiles/*/*_R2*_*.fa;do sed -i 's|^\(>.*\) |\1\\2 |' $f;done


## QC: count the target reads
wc -l targetReads_grep_plasmidClone/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $1, $4}' > targetReads_grep_plasmidClone/countsA
wc -l targetReads_grep_plasmidClone/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $5}' | awk -F'_' '{print $1,$2,$3,$4}' > targetReads_grep_plasmidClone/countsB
paste -d" " targetReads_grep_plasmidClone/countsA targetReads_grep_plasmidClone/countsB > targetReads_grep_plasmidClone/counts.list 
cat targetReads_grep_plasmidClone/counts.list | awk '{A[$3]+=$1}END{for(i in A)print i,A[i]}' > targetReads_grep_plasmidClone/count_per_sample.list

## merge target reads files
mkdir targetReads_grep_plasmidClone/mergedFiles/
while read sample count;do
  echo $sample;
  cat targetReads_grep_plasmidClone/singleFiles/*/${sample}_*.fa > targetReads_grep_plasmidClone/mergedFiles/${sample}.fa
done < targetReads_grep_plasmidClone/count_per_sample.list
#for sample in AVEAY003A AVEAY003B AVEAY001B AVEAY001A AVEAY15B;do
#cat targetReads_grep_plasmidClone/singleFiles/*/${sample}_*.fa > targetReads_grep_plasmidClone/mergedFiles/${sample}.fa
#done

## Multiple sequence alignment
#module load ClustalO/1.1.0
#cat polled.fa horned.fa target_reads.fa > all_seq.fa
#clustalo -i all_seq.fa -o all_seq_align.fa -t DNA -v

## sequence alignment to the whole genome
mkdir targetReads_grep_plasmidClone/mapped
module swap OpenMPI  OpenMPI/2.1.1
module load BWA/0.7.17
for f in targetReads_grep_plasmidClone/mergedFiles/*.fa;do
#for f in targetReads_grep_plasmidClone/mergedFiles/{AVEAY003A,AVEAY003B,AVEAY001B,AVEAY001A,AVEAY15B}.fa;do
  f2=targetReads_grep_plasmidClone/mapped/$(basename $f)
  bwa mem refGenome/BwaIndex/genome $f > ${f2%.fa}.genome.sam
done
cat targetReads_grep_plasmidClone/mapped/AVEAY15B.genome.sam |grep -v "^@" |awk '{print $2}' | sort | uniq -c 


#module swap GNU GNU/5.4.0-2.26
#module swap OpenMPI OpenMPI/1.10.3
#module load SAMtools/1.5
#for f in targetReads_grep/mergedFiles/*.fa;do
#  f2=targetReads_grep/mapped/$(basename $f)
#  samtools sort -T ${f2%.fa}_tempgen -O bam ${f2%.fa}.genome.sam > ${f2%.fa}.genome.sorted.bam
#  samtools index ${f2%.fa}.genome.sorted.bam
#done

## calculate the median kmer count normalized total kmer count and unique kmer count (instead of the read count 
# use ntcard to calc unique(u) and total(t) kmer count
# normalized as follow: 1. kmer count of grep reads and calc  median kmer count (m). 2. normalize m' = m/t * u/t * normlization factor 

## exclude non-perfectly matching reads (it would be great if we can allow up to 2 single base changes)  
# 1. map the reads to the bait
mkdir targetReads_grep_plasmidClone/bait_mapped
module swap OpenMPI  OpenMPI/2.1.1
module load BWA/0.7.17
bwa index targetReads_grep_plasmidClone/bait_unwrapped.fa
for f in targetReads_grep_plasmidClone/mergedFiles/*.fa;do
#for f in targetReads_grep_plasmidClone/mergedFiles/{AVEAY003A,AVEAY003B,AVEAY001B,AVEAY001A,AVEAY15B}.fa;do
  f2=targetReads_grep_plasmidClone/bait_mapped/$(basename $f)
  bwa mem targetReads_grep_plasmidClone/bait_unwrapped.fa $f > ${f2%.fa}.bait.sam
done
cat targetReads_grep_plasmidClone/bait_mapped/AVEAY{15B,003A,004A,003B,001B}.bait.sam |grep -v "^@" |awk '{print $2}' | sort | uniq -c  ## 575 0, 600 16, 18 2064, 4 4
cat targetReads_grep_plasmidClone/bait_mapped/AVEAY{007B,011A,012B,013A,14A,14B}.bait.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c ## 8 0, 5 16


mkdir -p targetReads_grep_plasmidClone/suspcious_reads && cd targetReads_grep_plasmidClone
for sample in AVEAY15B;do
  cat mapped/$sample.genome.sam |grep -v "^@" |awk '{if(($2==2048 || $2==2064) && ($5 >= 50))print $1}' > suspcious_reads/${sample}_genome.lst
  cat suspcious_reads/${sample}_genome.lst | grep -Fwf - mapped/$sample.genome.sam | cut -f1,2,3,4,5,6,7,8,9 > suspcious_reads/$sample.genome_only.sam
  cat suspcious_reads/$sample.genome_only.sam | awk '{if($2==0 || $2==16)print $0;}' > suspcious_reads/$sample.genome_only_1ry.sam
  cat suspcious_reads/$sample.genome_only.sam | awk '{if($2==2048 || $2==2064)print $0;}' > suspcious_reads/$sample.genome_only_sup.sam
  cat mapped/$sample.genome.sam |grep -v "^@" |awk '{if($4>=2428490 && $4<=2429901 && $6!="151M")print $0}'  > suspcious_reads/$sample.genome_only_del.sam

  cat bait_mapped/$sample.bait.sam |grep -v "^@" |awk '{if(($2==2048 || $2==2064) && ($5 >= 50))print $1}' > suspcious_reads/${sample}_bait.lst
  cat suspcious_reads/${sample}_bait.lst | grep -Fwf - bait_mapped/$sample.bait.sam | cut -f1,2,3,4,5,6,7,8,9 > suspcious_reads/$sample.bait_only.sam
  cat suspcious_reads/${sample}_genome.lst suspcious_reads/${sample}_bait.lst | sort | uniq > suspcious_reads/${sample}_both.lst
  cat suspcious_reads/${sample}_both.lst | grep -Fwf - mapped/$sample.genome.sam | cut -f1,2,3,4,5,6,7,8,9 > suspcious_reads/$sample.genome.sam
  cat suspcious_reads/${sample}_both.lst | grep -Fwf - bait_mapped/$sample.bait.sam | cut -f1,2,3,4,5,6,7,8,9 > suspcious_reads/$sample.bait.sam
done


###################################################

# 2. explore the mapping flag in all samples
cat targetReads_grep_plasmidClone/mapped/*.genome.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c  ##  18 0, 25 16, 1428 4
cat targetReads_grep_plasmidClone/mapped/AVEAY{15B,003A,004A,003B,001B}.genome.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c  ## 14 0, 22 16, 1143 4
cat targetReads_grep_plasmidClone/mapped/AVEAY{007B,011A,012B,013A,14A,14B}.genome.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c ## 2 16,11 4
# 3. characterize the mapping behiour of each read in bait and genome
mkdir -p targetReads_grep_plasmidClone/compare_map/edit
for f in targetReads_grep_plasmidClone/bait_mapped/AVEAY{15B,003A,004A,003B,001B}.bait.sam;do
  newf=targetReads_grep_plasmidClone/compare_map/edit/$(basename $f).info
  grep -v ^@ $f  | cut -f1,2,3,4,5,6,7,8,9  > $newf
  f2=targetReads_grep_plasmidClone/mapped/$(basename ${f%.bait.sam}).genome.sam
  newf2=targetReads_grep_plasmidClone/compare_map/edit/$(basename $f2).info
  grep -v ^@ $f2  | cut -f1,2,3,4,5,6,7,8,9  > $newf2
done

mkdir -p targetReads_grep_plasmidClone/compare_map/nonedit
for f in targetReads_grep_plasmidClone/bait_mapped/AVEAY{007B,011A,012B,013A,14A,14B}.bait.sam;do
  newf=targetReads_grep_plasmidClone/compare_map/nonedit/$(basename $f).info
  grep -v ^@ $f  | cut -f1,2,3,4,5,6,7,8,9  > $newf
  f2=targetReads_grep_plasmidClone/mapped/$(basename ${f%.bait.sam}).genome.sam
  newf2=targetReads_grep_plasmidClone/compare_map/nonedit/$(basename $f2).info
  grep -v ^@ $f2  | cut -f1,2,3,4,5,6,7,8,9  > $newf2
done


