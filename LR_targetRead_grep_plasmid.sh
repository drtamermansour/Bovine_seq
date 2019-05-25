work_dir=$(pwd)
mkdir -p ./targetReads_grep_plasmid
## data processing
# ....

## find target reads
cd $work_dir/targetReads_grep_plasmid ## The folder where you have the bait.fa files

## unwrap the bait sequences 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < bait.fa | tail -n +2 > bait_unwrapped.fa

## convert the bait to kmers
k=21
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
mkdir -p targetReads_grep_plasmid/singleFiles/{fastq_fail,fastq_pass}
while read name;do
  read=${name%%.fastq};
  echo $read
  seq="data/${read}.fastq"
  bait="targetReads_grep_plasmid/bait_Allkmers.txt"
  fq="targetReads_grep_plasmid/singleFiles/${read}_target_reads.fq"
  fa="targetReads_grep_plasmid/singleFiles/${read}_target_reads.fa"
  grep -B1 -A2 -F -f ${bait} ${seq} | grep -v "^\-\-" > ${fq}
  #cat ${fq} | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > ${fa}
  zgrep -B1 -F -f ${bait} ${seq} | grep -v "^\-\-" | sed 's/^@/>/' > ${fa}
done < <(ls data/fastq_*/*.fastq | sed 's|data/||')


## QC: count the target reads
wc -l targetReads_grep_plasmid/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $1, $4}' > targetReads_grep_plasmid/countsA
#wc -l targetReads_grep_plasmid/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $5}' | awk -F'_' '{print $1,$2,$3,$4}' > targetReads_grep_plasmid/countsB
wc -l targetReads_grep_plasmid/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $5}' | cut -d'_' -f1,2,3 > targetReads_grep_plasmid/countsB
paste -d" " targetReads_grep_plasmid/countsA targetReads_grep_plasmid/countsB > targetReads_grep_plasmid/counts.list 
cat targetReads_grep_plasmid/counts.list | awk '{A[$2]+=$1}END{for(i in A)print i,A[i]}' > targetReads_grep_plasmid/count_per_sample.list

## merge target reads files
mkdir targetReads_grep_plasmid/mergedFiles/
while read sample count;do
  echo $sample;
  cat targetReads_grep_plasmid/singleFiles/${sample}/*.fa > targetReads_grep_plasmid/mergedFiles/${sample}.fa
done < targetReads_grep_plasmid/count_per_sample.list
#############################################################################


