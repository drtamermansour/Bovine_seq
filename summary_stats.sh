# trimming report
echo "sample_ids" > reports/trim_ids
echo -e "Input_Read_Pairs\tBoth_Surviving\t%\tForward_Only_Surviving\t%\tReverse_Only_Surviving\t%\tDropped\t%" > reports/trim
for f in logs/snakejob.trimmomatic_pe.*.e*;do if [ "$(grep "Completed successfully" $f)" != "" ];then 
  grep "input: " $f | awk 'BEGIN{FS=OFS="/"}{print $6,$7}' | sed 's/_R2.*//' >> reports/trim_ids
  grep "Input Read Pairs: " $f | awk 'BEGIN{OFS="\t"}{print $4, $7, $8, $12, $13, $17, $18, $20, $21}';fi;done >> reports/trim
paste reports/trim_ids reports/trim > reports/trim_report
rm reports/trim_ids reports/trim

#coverage
> reports/cov_ids.txt
for f in qc/mappingQC/*.cov;do 
  echo $(basename $f) >> reports/cov_ids.txt 
  cat $f;done > reports/cov_stat.txt
paste reports/cov_ids.txt reports/cov_stat.txt | awk '{print $1,$4}' > reports/cov.txt
rm reports/cov_ids.txt reports/cov_stat.txt

#flagstat
echo "mapping_ids" > reports/mapping_ids;
cat qc/mappingQC/AVEAY0010A.stat | sed 's/ (.*)//' | cut -d" " --complement -f1,2,3 | sed 's/ /_/g' | tr "\n" "\t" > reports/mapping_report
echo "" >> reports/mapping_report
for f in qc/mappingQC/*.stat;do 
  echo $(basename $f) >> reports/mapping_ids; 
  awk '{print $1}' $f | tr "\n" "\t";echo "";done >> reports/mapping_report
paste reports/mapping_ids reports/mapping_report > reports/mapping_report_ids
echo -e "no_of_reads(read1+read2)\tMapped\tSecondary_alignments\tBoth_mapped\tProperly_paired\tmate_mapped_to_different_chr\tsingletons" > reports/mapping_percentages
tail -n+2 reports/mapping_report | awk 'BEGIN{FS=OFS="\t"}{print $6, $5/$1, $2/$1, $10/$6, $9/$6, $12/$6, $11/$6}' >> reports/mapping_percentages
paste reports/mapping_ids reports/mapping_percentages > reports/mapping_percentages_ids
rm reports/mapping_ids reports/mapping_percentages reports/mapping_report

#dedub metrics
echo "mapping_ids" > reports/dedup_ids;
tail -n6 data/dedup/AVEAY001A.txt | head -n1 > reports/dedup_summary
for f in data/dedup/*.txt;do 
  for i in {1..3};do echo $(basename $f);done >> reports/dedup_ids  
  tail -n5 $f | head -n3;done >> reports/dedup_summary 
paste reports/dedup_ids reports/dedup_summary > reports/dedup_summary_ids
rm reports/dedup_ids reports/dedup_summary

#variant calling stats
#https://software.broadinstitute.org/gatk/documentation/article?id=11069
input="vc/hapCaller_raw.vcf"
mkdir -p vc/rawQC
for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AC" "AF" "AN" "InbreedingCoeff";do
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > vc/rawQC/$filter.report
done

for report in vc/rawQC/*.report;do
  Rscript -e 'args=(commandArgs(TRUE)); input=args[1]; data=read.table(input);'\
'library(ggplot2); outputPDF=paste(input,"pdf",sep=".");'\
'pdf(outputPDF);ggplot(data, aes(data$V2)) + geom_density(alpha = 0.2);dev.off();' $report
done

awk '{A[$2]++}END{for(i in A)print i,A[i]}' vc/rawQC/AN.report | sort -n > vc/rawQC/AN.report.summary 
