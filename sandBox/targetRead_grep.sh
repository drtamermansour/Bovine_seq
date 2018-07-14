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
done < bait_unwrapped.fa | tr "atcg" "ATCG" | sort | uniq > bait_kmers.txt

## find reverse complement kmers
cat bait_kmers.txt | rev | tr "ATGC" "TACG" > bait_rckmers.txt

## make one list of all possible kmers
cat bait_kmers.txt bait_rckmers.txt | sort | uniq > bait_Allkmers.txt

## find matching reads (in Snakemake)
cd $work_dir
mkdir targetReads_grep/singleFiles
data="data/trimmed"
output="targetReads_grep/singleFiles"
for f in $data/*/*.fastq.gz;do
 batch=$(basename $(dirname $f))
 infile=$(basename $f) 
 outbase=${infile%.fastq.gz}_target_reads
 zgrep -B1 -A2 -F -f targetReads_grep/bait_Allkmers.txt $f | grep -v "^\-\-" > $output/$batch/$outbase.fq
 cat $output/$batch/$outbase.fq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $output/$batch/$outbase.fa
done

## recognize the ends
for f in $output/*/*_R1_*.fa;do sed -i 's|^\(>.*\) |\1\\1 |' $f;done 
for f in $output/*/*_R2_*.fa;do sed -i 's|^\(>.*\) |\1\\2 |' $f;done

## QC: count the target reads
wc -l targetReads_grep/singleFiles/1*/*.fa | head -n400 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $1, $3}' > targetReads_grep/countsA
wc -l targetReads_grep/singleFiles/1*/*.fa | head -n400 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $4}' | awk -F'_' '{print $1,$3,$4}' > targetReads_grep/countsB
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
module load bwa/0.7.7.r441
module load SAMTools/1.5
for f in targetReads_grep/mergedFiles/*.fa;do
  bwa mem refGenome/BwaIndex/genome $f > ${f%.fa}.genome.sam
  samtools sort -T ${f%.fa}_tempgen -O bam ${f%.fa}.genome.sam > ${f%.fa}.genome.sorted.bam
  samtools index ${f%.fa}.genome.sorted.bam
done

## get the reads back from the sam to have the same strand
f="targetReads_grep/mergedFiles/AVEAY001A.genome.sam" ## f="targetReads_grep/mergedFiles/AVEAY0010B.genome.sam"
grep -v "@" $f | awk 'BEGIN{FS="\t";OFS="\n"}{print ">"$1,$10}' > $f.fa
## assemble by extension 
/opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $f.fa -w 1 -m 16 -o 1;

## error trimming and assembly by extention 
cd targetReads_grep/mergedFiles/
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
for f in targetReads_grep/mergedFiles/*.fa;do
  bwa mem targetReads_grep/bait.fa $f > ${f%.fa}.sam
  samtools sort -T ${f%.fa}_temp -O bam ${f%.fa}.sam > ${f%.fa}.sorted.bam
  samtools index ${f%.fa}.sorted.bam
done
#samtools tview target_reads_sorted.bam polled.fa

## get the reads back from the sam to have the same strand
f="targetReads_grep/mergedFiles/AVEAY001A.sam" ## "targetReads_grep/mergedFiles/AVEAY010B.sam"
grep -v "@" $f | awk 'BEGIN{FS="\t";OFS="\n"}{if($12=="NM:i:0" && $13=="MD:Z:151")print ">"$1,$10}' > $f.per.fa
grep -v "@" $f | awk 'BEGIN{FS="\t";OFS="\n"}{if($12!="NM:i:0" || $13!="MD:Z:151")print ">"$1,$10}' > $f.mis.fa
## assemble by extension 
/opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $f.per.fa -w 1 -m 16 -o 1;
/opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $f.mis.fa -w 2 -m 16 -o 1;


##############################
## create new_bait.txt for the read with unknown insertion then
for f in data/trimmed/*/AVEAY001A*.gz;do zgrep -B1 -A2 -F -f targetReads_grep/new_bait.txt $f | grep -v "^\-\-" >> targetReads_grep/new_bait.fq;done 
cat targetReads_grep/new_bait.fq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > targetReads_grep/new_bait.fa


