work_dir=$(pwd)
mkdir -p ./{data,targetReads_grep}
## data processing
# ....

## find target reads
cd $work_dir/targetReads_grep ## The folder where you have the bait.fa, polled.fa, and horned.fa files

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

## find matching reads (in Snakemake)
#cd $work_dir
#mkdir targetReads_grep/singleFiles
#data="data/trimmed"
#output="targetReads_grep/singleFiles"
#for f in $data/*/*.fastq.gz;do
# batch=$(basename $(dirname $f))
# infile=$(basename $f) 
# outbase=${infile%.fastq.gz}_target_reads
# zgrep -B1 -A2 -F -f targetReads_grep/bait_Allkmers.txt $f | grep -v "^\-\-" > $output/$batch/$outbase.fq
# cat $output/$batch/$outbase.fq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $output/$batch/$outbase.fa
#done

## recognize the ends
for f in targetReads_grep/singleFiles/*/*_R1*_*.fa;do sed -i 's|^\(>.*\) |\1\\1 |' $f;done
for f in targetReads_grep/singleFiles/*/*_R2*_*.fa;do sed -i 's|^\(>.*\) |\1\\2 |' $f;done
#for f in targetReads_grep/singleFiles/{1_4_18-1,1_4_18-2,12_19_17}/*_R1*_*.fa;do sed -i 's|^\(>.*\) |\1\\1 |' $f;done
#for f in targetReads_grep/singleFiles/{1_4_18-1,1_4_18-2,12_19_17}/*_R2*_*.fa;do sed -i 's|^\(>.*\) |\1\\2 |' $f;done
#for f in targetReads_grep/singleFiles/{5_2_18,6_13_18,8_2_18,10_9_18,10_17_18}/*_R1*_*.fa;do sed -i 's|^\(>.*\) |\1\\1 |' $f;done 
#for f in targetReads_grep/singleFiles/{5_2_18,6_13_18,8_2_18,10_9_18,10_17_18}/*_R2*_*.fa;do sed -i 's|^\(>.*\) |\1\\2 |' $f;done

## QC: count the target reads
wc -l targetReads_grep/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $1, $4}' > targetReads_grep/countsA
wc -l targetReads_grep/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $5}' | awk -F'_' '{print $1,$2,$3,$4}' > targetReads_grep/countsB
paste -d" " targetReads_grep/countsA targetReads_grep/countsB > targetReads_grep/counts.list 
cat targetReads_grep/counts.list | awk '{A[$3]+=$1}END{for(i in A)print i,A[i]}' > targetReads_grep/count_per_sample.list

## merge target reads files
mkdir targetReads_grep/mergedFiles/
while read sample count;do
  echo $sample;
  cat targetReads_grep/singleFiles/*/${sample}_*.fa > targetReads_grep/mergedFiles/${sample}.fa
done < targetReads_grep/count_per_sample.list

## Multiple sequence alignment
#module load ClustalO/1.1.0
#cat polled.fa horned.fa target_reads.fa > all_seq.fa
#clustalo -i all_seq.fa -o all_seq_align.fa -t DNA -v

## sequence alignment to the whole genome
mkdir targetReads_grep/mapped
module swap OpenMPI  OpenMPI/2.1.1
module load BWA/0.7.17
for f in targetReads_grep/mergedFiles/*.fa;do
  f2=targetReads_grep/mapped/$(basename $f)
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
# 1. find the coordinates of the 212bp in the new assembly
module swap OpenMPI  OpenMPI/2.1.1
module load BWA/0.7.17
bwa mem refGenome/BwaIndex/genome targetReads_grep/bait.fa > bait.genome.sam #CM008168.2:2429109- 2429320  
# 2. explore the mapping flag in all samples
cat targetReads_grep/mapped/*.genome.sam | grep -v "^@" | awk '{print $2}' | sort | uniq -c
# 3. exclude all reads normally mapping to the 212bp locus (add 150 bp before and after to consider the hanging reads
mkdir -p targetReads_grep/abnReads_temp
for f in targetReads_grep/mapped/*.genome.sam;do
  f2=targetReads_grep/abnReads_temp/$(basename $f)
  cat $f | grep -v "^@" | awk '{if($4< 2428959 || $4> 2429320 || ($2!=0 && $2!=16)) print $0}' > $f2
done
find targetReads_grep/abnReads_temp -size  0 -print0 |xargs -0 rm --
mkdir -p targetReads_grep/abnReads
for f in targetReads_grep/abnReads_temp/*.genome.sam;do
  f2=targetReads_grep/mapped/$(basename $f)
  f3=targetReads_grep/abnReads/$(basename $f)
  grep "^@" $f2 > $f3
  cat $f >> $f3
done  
## get the reads back
for f in targetReads_grep/abnReads/*.genome.sam;do
  #grep -v "@" $f | awk 'BEGIN{FS="\t";OFS="\n"}{print ">"$1,$10}' > $f.fa
  f2=targetReads_grep/mergedFiles/$(basename $f .genome.sam).fa
  grep -v "@" $f | awk '{print $1}' | sort | uniq | grep -Fwf - -A1 $f2 | grep -v "^--$" > $f.fa
  wc -l $f.fa
done

## sequence alignment to the target insertion
cat external/*/targetReads_grep/abnReads/*.genome.sam.fa targetReads_grep/abnReads/*.genome.sam.fa > abnReads.fa
module swap OpenMPI  OpenMPI/2.1.1
module load BWA/0.7.17
bwa index targetReads_grep/insertion.fa
bwa mem targetReads_grep/insertion.fa abnReads.fa > abnReads.insertion.sam
module swap GNU GNU/5.4.0-2.26
module swap OpenMPI OpenMPI/1.10.3
module load SAMtools/1.5
samtools sort -T abnReads_tempgen -O bam abnReads.insertion.sam > abnReads.insertion.sorted.bam
samtools index abnReads.insertion.sorted.bam


## 4. select all reads normally mapping to the 212bp locus (add 150 bp before and after to consider the hanging reads
mkdir -p targetReads_grep/normReads
for f in targetReads_grep/mapped/*.genome.sam;do
  f2=targetReads_grep/normReads/$(basename $f .genome.sam)
  f3=targetReads_grep/abnReads_temp/$(basename $f)
  f4=targetReads_grep/normReads/$(basename $f).temp
  f5=targetReads_grep/mergedFiles/$(basename $f .genome.sam).fa
  if [ -f $f3 ];then awk '{print $1}' $f3 > $f4;else echo "" > $f4;fi
  grep -v "@" $f | awk '{print $1}' | sort | uniq > $f4.2
  grep -vFwf $f4 $f4.2 | grep -Fwf - -A1 $f5 | grep -v "^--$" > $f2.fa
  wc -l $f2.fa
  rm ${f4}*
done

cat external/*/targetReads_grep/normReads/SRR3290{535,615}.fa targetReads_grep/normReads/AVEAY{15B,001A,001B,002A,002B,003A,003B,004A,004B}.fa > normReads_pc.fa
bwa mem refGenome/BwaIndex/genome normReads_pc.fa > normReads_pc.genome.sam
cat external/*/targetReads_grep/normReads/SRR329063{1,2}.fa targetReads_grep/normReads/AVEAY{1[456]A,14B,01[123][AB],00[56789][AB],0010[AB]}.fa > normReads_p.fa
bwa mem refGenome/BwaIndex/genome normReads_p.fa > normReads_p.genome.sam
module swap GNU GNU/5.4.0-2.26
module swap OpenMPI OpenMPI/1.10.3
module load SAMtools/1.5
samtools sort -T normReads_pc_tempgen -O bam normReads_pc.genome.sam > normReads_pc.genome.sorted.bam
samtools index normReads_pc.genome.sorted.bam
samtools sort -T normReads_p_tempgen -O bam normReads_p.genome.sam > normReads_p.genome.sorted.bam
samtools index normReads_p.genome.sorted.bam



# assembl bynextension 
#source activate workEnv1
#for f in targetReads_grep/abnReads/*.genome.sam.fa;do
#  SSAKE -f $f -w 1 -m 16 -o 1;
#done 
##########################
## TEMP code
## get the reads back from the sam to have the same strand
mkdir targetReads_grep/mergedFiles_oneStrand
f="targetReads_grep/mergedFiles/AVEAY001A.genome.sam" ## f="targetReads_grep/mergedFiles/AVEAY0010B.genome.sam"
f2=targetReads_grep/mergedFiles_oneStrand/$(basename $f)
grep -v "@" $f | awk 'BEGIN{FS="\t";OFS="\n"}{print ">"$1,$10}' > $f2.fa
## assemble by extension 
/opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $f2.fa -w 1 -m 16 -o 1;

## error trimming and assembly by extention 
cd targetReads_grep/mergedFiles_oneStrand/
module swap GNU GNU/4.8.2
module load khmer/2.0
f="AVEAY001A.genome.sam" ## f="AVEAY0010B.genome.sam"
load-into-counting.py -k 20 -x 5e7 $f.countgraph $f.fa
filter-abund.py -C 2 $f.countgraph $f.fa
#processed 59 / wrote 54 / removed 5                   # processed 26 / wrote 24 / removed 2
#processed 8602 bp / wrote 7468 bp / removed 1134 bp   # processed 3926 bp / wrote 3033 bp / removed 893 bp
#discarded 13.2%                                       # discarded 22.7%
/opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $f.fa.abundfilt -w 1 -m 16 -o 1;


## classify reads into perfect matching and mismatching
#f="targetReads_grep/mergedFiles/AVEAY001A.genome.sam"
#grep -v "@" $f | awk '{if($12=="NM:i:0" && $13=="MD:Z:151")print $0}' > $f.Pmatch
#grep -v "@" $f | awk '{if($12!="NM:i:0" || $13!="MD:Z:151")print $0}' > $f.Mmatch

##############################
## sequence alignment to the target reference
module load bwa/0.7.7.r441
module load SAMTools/1.5
#bwa index targetReads_grep/polled.fa
bwa index targetReads_grep/bait.fa
mkdir targetReads_grep/bait_mapped
for f in targetReads_grep/mergedFiles/*.fa;do
  f2=targetReads_grep/bait_mapped/$(basename $f)
  bwa mem targetReads_grep/bait.fa $f > ${f2%.fa}.sam
  samtools sort -T ${f2%.fa}_temp -O bam ${f2%.fa}.sam > ${f2%.fa}.sorted.bam
  samtools index ${f2%.fa}.sorted.bam
done
#samtools tview target_reads_sorted.bam polled.fa

## get the reads back from the sam to have the same strand
mkdir targetReads_grep/mergedFiles_baitStrand
f="targetReads_grep/bait_mapped/AVEAY001A.sam" ## "targetReads_grep/bait_mapped/AVEAY010B.sam"
f2=targetReads_grep/mergedFiles_baitStrand/$(basename $f)
grep -v "@" $f | awk 'BEGIN{FS="\t";OFS="\n"}{if($12=="NM:i:0" && $13=="MD:Z:151")print ">"$1,$10}' > $f2.per.fa
grep -v "@" $f | awk 'BEGIN{FS="\t";OFS="\n"}{if($12!="NM:i:0" || $13!="MD:Z:151")print ">"$1,$10}' > $f2.mis.fa
## assemble by extension 
/opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $f.per.fa -w 1 -m 16 -o 1;
/opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $f.mis.fa -w 2 -m 16 -o 1;


##############################
## create new_bait.txt for the read with unknown insertion then
for f in data/trimmed/*/AVEAY001A*.gz;do zgrep -B1 -A2 -F -f targetReads_grep/new_bait.txt $f | grep -v "^\-\-" >> targetReads_grep/new_bait.fq;done 
cat targetReads_grep/new_bait.fq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > targetReads_grep/new_bait.fa


