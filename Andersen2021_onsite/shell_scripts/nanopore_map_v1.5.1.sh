#!/bin/bash 
module load Minimap2/2.17-foss-2018a
module load SAMtools/1.10-foss-2018a

###################################################################################################
#
#  Version 1.5.1
#
#  This workflow trims Nanopore readnames by removing runID, read# and read length.
#  The rest of the information (readID, read start time and barcode) are retained in the header.
#  The resulting fastq file is mapped to a reference database and further processed to retain
#  only the readID, read start time, barcode, the matched OTU and the mapping identity
#
#  Martin Andersen, 2019
#
###################################################################################################

READS=seqdata/nanopore/gup315_called_demulti/np_all_gup315_trim_filt;
mapper=minimap2;
TODAY=$(date '+%Y-%m-%d')
DATABASE=databases/midas37/midas37_99clust.fa;

clear
echo ""
if [ $mapper = "minimap2" ] || [ $mapper = "usearch10" ] ; then
	echo "Running nanopore read mapping version 1.5.1, using $mapper"
else
	echo "The only supported mappers are usearch10 and minimap2"
	echo ""
	exit 1
fi
date

echo ""
echo "Trimming read names, removing runID, read# and length"
ENDING=".fastq";
READSEND=$READS$ENDING;
TRIM="_nametrim.fastq";
READIDTRIM=$READS$TRIM;
sed 's/\(runid.*start_\)//' < $READSEND | tr -d '[:blank:]' > seqdata/nanopore/gup315_called_demulti/np_all_gup315_trim_filt_nametrim.fastq
date

echo ""
echo "Shortened read names. Starting mapping to database"
sam=".sam";
MAPPINGNAME=$MAPNAME$sam;
if [ $mapper = "minimap2" ]
then
	minimap2 -ax map-ont -t 40 --secondary=no --MD $DATABASE seqdata/nanopore/gup315_called_demulti/np_all_gup315_trim_filt_nametrim.fastq > mappings/minimap-npv13-lib4-gup315-midas37.sam
elif [ $mapper = "usearch10" ]
then
	usearch10 -usearch_local $READS -db $DATABASE -strand both -id 0.8 -top_hit_only -maxaccepts 1 -samout $MAPPINGNAME -threads 20
else
	echo "The only supported mapping tools are usearch10 or minimap2"
	exit 1
fi
date

if [ $mapper = "minimap2" ]
then
	echo ""
	echo "Mapped reads. Removing supplementary mappings based on SAM bit flags"
	samtools view -F 256 -F 4 -F 2048 mappings/minimap-npv13-lib4-gup315-midas37.sam -o mappings/minimap-npv13-lib4-gup315-midas37_nodupes.sam
	date
else
	:
fi

sed '/^@/ d' mappings/minimap-npv13-lib4-gup315-midas37_nodupes.sam | \
	awk '{
    for(i=1;i<=NF;i++){
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); mm = $i}
    }
    split($6, count, /[^0-9]+/);
    split($6, type, /[^A-Z]*/);
    for(i=1; i<= length(count)-1; i++){
      if(type[i + 1] ~ /[DIM]/){aln+=count[i]};
    }
    print $1, $2, $3, length($10), aln, (aln - mm)/aln, $12, $14, $20
    aln=0;
  }' > np_lib4_idmapped.txt

sed 's/barcode/\tbarcode/' np_lib4_idmapped.txt | sed 's/time=/\t/' > mappings/${TODAY}_minimap-npv13-lib4-gup315-midas37-99clust.txt