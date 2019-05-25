work_dir=$(pwd)
mkdir -p ./targetReads_grep_plasmid
## data processing
# ....

## find target reads
cd $work_dir/targetReads_grep_plasmid ## The folder where you have the bait.fa files

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
mkdir -p targetReads_grep_plasmid/singleFiles/{10_17_18,10_9_18,12_19_17,1_4_18-1,1_4_18-2,5_2_18,6_13_18,8_2_18}
while read name;do
  read=${name%%.fastq.gz};
  echo $read
  seq="data/trimmed/${read}.fastq.gz"
  bait="targetReads_grep_plasmid/bait_Allkmers.txt"
  fq="targetReads_grep_plasmid/singleFiles/${read}_target_reads.fq"
  fa="targetReads_grep_plasmid/singleFiles/${read}_target_reads.fa"
  #zgrep -B1 -A2 -F -f ${bait} ${seq} | grep -v "^\-\-" > ${fq}
  #cat ${fq} | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > ${fa}
  zgrep -B1 -F -f ${bait} ${seq} | grep -v "^\-\-" | sed 's/^@/>/' > ${fa}
done < <(ls data/fastq/*/*.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/*/AVEAY003B_*.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/*/AVEAY003A_*.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/10_17_18/AVEAY15B_S14_L006_R2_001.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/*/AVEAY001A_*.fastq.gz | sed 's|data/fastq/||')
#done < <(ls data/fastq/1_4_18-1/AVEAY001B_S2_L007_R2_001.fastq.gz | sed 's|data/fastq/||')

## recognize the ends
for f in targetReads_grep_plasmid/singleFiles/*/*_R1*_*.fa;do sed -i 's|^\(>.*\) |\1\\1 |' $f;done
for f in targetReads_grep_plasmid/singleFiles/*/*_R2*_*.fa;do sed -i 's|^\(>.*\) |\1\\2 |' $f;done


## QC: count the target reads
wc -l targetReads_grep_plasmid/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $1, $4}' > targetReads_grep_plasmid/countsA
wc -l targetReads_grep_plasmid/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $5}' | awk -F'_' '{print $1,$2,$3,$4}' > targetReads_grep_plasmid/countsB
paste -d" " targetReads_grep_plasmid/countsA targetReads_grep_plasmid/countsB > targetReads_grep_plasmid/counts.list 
cat targetReads_grep_plasmid/counts.list | awk '{A[$3]+=$1}END{for(i in A)print i,A[i]}' > targetReads_grep_plasmid/count_per_sample.list

## merge target reads files
mkdir targetReads_grep_plasmid/mergedFiles/
while read sample count;do
  echo $sample;
  cat targetReads_grep_plasmid/singleFiles/*/${sample}_*.fa > targetReads_grep_plasmid/mergedFiles/${sample}.fa
done < targetReads_grep_plasmid/count_per_sample.list
#for sample in AVEAY003A AVEAY003B AVEAY001B AVEAY001A AVEAY15B;do
#cat targetReads_grep_plasmid/singleFiles/*/${sample}_*.fa > targetReads_grep_plasmid/mergedFiles/${sample}.fa
#done

## Multiple sequence alignment
#module load ClustalO/1.1.0
#cat polled.fa horned.fa target_reads.fa > all_seq.fa
#clustalo -i all_seq.fa -o all_seq_align.fa -t DNA -v

## sequence alignment to the whole genome
mkdir targetReads_grep_plasmid/mapped
module swap OpenMPI  OpenMPI/2.1.1
module load BWA/0.7.17
for f in targetReads_grep_plasmid/mergedFiles/*.fa;do
#for f in targetReads_grep_plasmid/mergedFiles/{AVEAY003A,AVEAY003B,AVEAY001B,AVEAY001A,AVEAY15B}.fa;do
  f2=targetReads_grep_plasmid/mapped/$(basename $f)
  bwa mem refGenome/BwaIndex/genome $f > ${f2%.fa}.genome.sam
done

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
mkdir targetReads_grep_plasmid/bait_mapped
module swap OpenMPI  OpenMPI/2.1.1
module load BWA/0.7.17
bwa index targetReads_grep_plasmid/bait.fa
for f in targetReads_grep_plasmid/mergedFiles/*.fa;do
#for f in targetReads_grep_plasmid/mergedFiles/{AVEAY003A,AVEAY003B,AVEAY001B,AVEAY001A,AVEAY15B}.fa;do
  f2=targetReads_grep_plasmid/bait_mapped/$(basename $f)
  bwa mem targetReads_grep_plasmid/bait.fa $f > ${f2%.fa}.bait.sam
done
cat targetReads_grep_plasmid/bait_mapped/AVEAY{15B,003A,004A,003B,001B}.bait.sam |grep -v "^@" |awk '{print $2}' | sort | uniq -c  ## 575 0, 600 16, 18 2064, 4 4
cat targetReads_grep_plasmid/bait_mapped/AVEAY{007B,011A,012B,013A,14A,14B}.bait.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c ## 8 0, 5 16
# 2. explore the mapping flag in all samples
cat targetReads_grep_plasmid/mapped/*.genome.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c  ##  18 0, 25 16, 1428 4
cat targetReads_grep_plasmid/mapped/AVEAY{15B,003A,004A,003B,001B}.genome.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c  ## 14 0, 22 16, 1143 4
cat targetReads_grep_plasmid/mapped/AVEAY{007B,011A,012B,013A,14A,14B}.genome.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c ## 2 16,11 4
# 3. characterize the mapping behaviour of each read in bait and genome
mkdir -p targetReads_grep_plasmid/compare_map/edit
for f in targetReads_grep_plasmid/bait_mapped/AVEAY{15B,003A,004A,003B,001B}.bait.sam;do
  newf=targetReads_grep_plasmid/compare_map/edit/$(basename $f).info
  grep -v ^@ $f  | cut -f1,2,3,4,5,6,7,8,9  > $newf
  f2=targetReads_grep_plasmid/mapped/$(basename ${f%.bait.sam}).genome.sam
  newf2=targetReads_grep_plasmid/compare_map/edit/$(basename $f2).info
  grep -v ^@ $f2  | cut -f1,2,3,4,5,6,7,8,9  > $newf2
done

mkdir -p targetReads_grep_plasmid/compare_map/nonedit
for f in targetReads_grep_plasmid/bait_mapped/AVEAY{007B,011A,012B,013A,14A,14B}.bait.sam;do
  newf=targetReads_grep_plasmid/compare_map/nonedit/$(basename $f).info
  grep -v ^@ $f  | cut -f1,2,3,4,5,6,7,8,9  > $newf
  f2=targetReads_grep_plasmid/mapped/$(basename ${f%.bait.sam}).genome.sam
  newf2=targetReads_grep_plasmid/compare_map/nonedit/$(basename $f2).info
  grep -v ^@ $f2  | cut -f1,2,3,4,5,6,7,8,9  > $newf2
done


## collect the the suspcious reads in a summary alignment file
echo "#AVEAY003B" > suspicious_reads.sam
echo "K00188:49:HN23WBBXX:7:2120:5984:35637\1" | grep -Fwf - targetReads_grep_plasmid/bait_mapped/AVEAY003B.bait.sam >> suspicious_reads.sam
echo "K00188:49:HN23WBBXX:7:2120:5984:35637\1" | grep -Fwf - targetReads_grep_plasmid/mapped/AVEAY003B.genome.sam >> suspicious_reads.sam
echo "#AVEAY004A" >> suspicious_reads.sam
echo "K00188:39:HMYWJBBXX:3:2210:25591:19144\1" | grep -Fwf - targetReads_grep_plasmid/bait_mapped/AVEAY004A.bait.sam >> suspicious_reads.sam
echo "K00188:39:HMYWJBBXX:3:2210:25591:19144\1" | grep -Fwf - targetReads_grep_plasmid/mapped/AVEAY004A.genome.sam >> suspicious_reads.sam
echo "#AVEAY007B" >> suspicious_reads.sam
echo "K00188:49:HN23WBBXX:7:1101:26768:10915\1" | grep -Fwf - targetReads_grep_plasmid/bait_mapped/AVEAY007B.bait.sam >> suspicious_reads.sam
echo "K00188:49:HN23WBBXX:7:1101:26768:10915\1" | grep -Fwf - targetReads_grep_plasmid/mapped/AVEAY007B.genome.sam >> suspicious_reads.sam
echo "#AVEAY012B" >> suspicious_reads.sam
echo "K00364:102:HWV25BBXX:6:2213:26058:16928\2" | grep -Fwf - targetReads_grep_plasmid/bait_mapped/AVEAY012B.bait.sam >> suspicious_reads.sam
echo "K00364:102:HWV25BBXX:6:2213:26058:16928\2" | grep -Fwf - targetReads_grep_plasmid/mapped/AVEAY012B.genome.sam >> suspicious_reads.sam

