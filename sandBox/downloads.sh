## download the data of the 1st 3 patches from the hard disk
mkdir Bovine_seq/data/fastq
cd Bovine_seq/data/fastq

## download the data from the 4th patch (ftp server
cd Bovine_seq/data/fastq/5_2_18/
for f in *.gz;do
newf=$(echo $f | sed 's/_/0/')
echo $newf;
mv $f $newf;
done

## ftp download
wget --no-verbose --no-parent --recursive --level=1 --no-directories --user=gslftp --password=gsl23ftp! ftp://gslserver.qb3.berkeley.edu/181004_150PE_HS4K2A/L235678/Young/AVEAY*

cd Bovine_seq/data/fastq/10_9_18/
for f in *.gz;do
newf=$(echo $f | sed 's/_/0/')
echo $newf;
mv $f $newf;
done


