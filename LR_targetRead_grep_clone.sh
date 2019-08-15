cd BovineSeq_LR
work_dir=$(pwd)
mkdir -p ./targetReads_grep_clone
## data processing
# ....

## find target reads
cd $work_dir/targetReads_grep_clone ## The folder where you have the bait.fa files

## unwrap the bait sequences 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < bait.fa | tail -n +2 > bait_unwrapped.fa

## convert the bait to kmers
k=31
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
mkdir -p targetReads_grep_clone/singleFiles/{fastq_fail,fastq_pass}
while read name;do
  read=${name%%.fastq};
  echo $read
  seq="data/${read}.fastq"
  bait="targetReads_grep_clone/bait_Allkmers.txt"
  fq="targetReads_grep_clone/singleFiles/${read}_target_reads.fq"
  fa="targetReads_grep_clone/singleFiles/${read}_target_reads.fa"
  grep -B1 -A2 --no-group-separator -F -f ${bait} ${seq} > ${fq}
  #cat ${fq} | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > ${fa}
  grep -B1 --no-group-separator -F -f ${bait} ${seq} | sed 's/^@/>/' > ${fa}
done < <(ls data/fastq_*/*.fastq | sed 's|data/||')


## QC: count the target reads
wc -l targetReads_grep_clone/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $1, $4}' > targetReads_grep_clone/countsA
#wc -l targetReads_grep_clone/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $5}' | awk -F'_' '{print $1,$2,$3,$4}' > targetReads_grep_clone/countsB
wc -l targetReads_grep_clone/singleFiles/*/*.fa | head -n -1 | sed -e 's/^[[:space:]]*//' | awk -F'[ /]' '{print $5}' | cut -d'_' -f1,2,3 > targetReads_grep_clone/countsB
paste -d" " targetReads_grep_clone/countsA targetReads_grep_clone/countsB > targetReads_grep_clone/counts.list 
cat targetReads_grep_clone/counts.list | awk '{A[$2]+=$1}END{for(i in A)print i,A[i]}' > targetReads_grep_clone/count_per_sample.list

## merge target reads files
mkdir targetReads_grep_clone/mergedFiles/
while read sample count;do
  echo $sample;
  cat targetReads_grep_clone/singleFiles/${sample}/*.fa > targetReads_grep_clone/mergedFiles/${sample}.fa
  cat targetReads_grep_clone/singleFiles/${sample}/*.fq > targetReads_grep_clone/mergedFiles/${sample}.fq
done < targetReads_grep_clone/count_per_sample.list


