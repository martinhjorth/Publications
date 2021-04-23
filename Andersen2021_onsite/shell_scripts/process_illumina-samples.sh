###################################################################################################
#
#  Version 4.3
#
#  This workflow script generates OTU tables from raw V13 16S amplicon data.
#
#  It is currently only supported for internal use on Aalborg University,
#  but feel free to use any parts that you find usefull.
#
#  Mads Albertsen, 2015
#
###################################################################################################

clear
echo ""
echo "Running: 16S V13 workflow version 4.3"
date

echo ""
echo "Finding your samples and copying them to the current folder"
while read samples
do
a="_";
NAME=$samples$a;
find /space/sequences/Illumina/MiSeq/ -name $NAME*R1* 2>/dev/null -exec gzip -cd {} \; | head -n 1  | sed 's/\@/>/' >> id.txt
find /space/sequences/Illumina/MiSeq/ -name $NAME*R1* 2>/dev/null -exec gzip -cd {} \; >> forward.fastq
find /space/sequences/Illumina/MiSeq/ -name $NAME*R2* 2>/dev/null -exec gzip -cd {} \; >> reverse.fastq
done < samples

paste -d "\t" id.txt samples > sampleid.txt
date

echo ""
echo "Removing bad quality reads"
trimmomatic PE -threads 20 forward.fastq reverse.fastq forward_qs.fastq.gz s1.fastq.gz reverse_qs.fastq.gz s2.fastq.gz SLIDINGWINDOW:5:3 MINLEN:275 2> temp.log
date

echo ""
echo "Merging reads"
flash -m 10 -M 200 forward_qs.fastq.gz reverse_qs.fastq.gz -o temp > flash.log
trim.fastq.length.pl -i temp.extendedFrags.fastq -o merged_l.fastq -m 425 -x 525 > temp.log
date

sed 's/ /Zbarcode=/' merged_l.fastq > il_lib5_libvar.fastq
