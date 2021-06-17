#!/bin/bash 
module load Minimap2/2.17-foss-2018a
module load SAMtools/1.10-foss-2018a

READS=seqdata/illumina/lib4/il_midasSTD-lib4.fastq;
TODAY=$(date '+%Y-%m-%d')

DB=db/midas37/midas37_99clust.fa;

module load Minimap2/2.17-foss-2018a
minimap2 -ax sr -t 40 --secondary=no --MD $DB $READS > output/mappings/minimap_il-lib4_midas37_99clust.sam

echo ""
echo "Mapped reads. Removing supplementary mappings based on SAM bit flags"
samtools view -F 256 -F 4 -F 2048 output/mappings/minimap_il-lib4_midas37_99clust.sam -o output/mappings/minimap_il-lib4_midas37_99clust-nodupes.sam
date

sed '/^@/ d' output/mappings/minimap_il-lib4_midas37_99clust-nodupes.sam | \
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
  }' > idmapped_ilv13_clust.txt

awk -F"\t" -v OFS='\t' '{print $1, $3}' idmapped_ilv13_clust.txt | sed 's/Zbarcode=/\t/' >  output/mappings/${TODAY}_minimap_il-lib4_midas37_99clust_processed.txt
