READS=seqdata/illumina/il_midasSTD.fastq;
DB=db/midas37/midas37_single.fa.fa;

# start mapping
module load Minimap2/2.17-foss-2018a
minimap2 -ax sr -t 20 --secondary=no --MD $DB $READS > map_db_test/mappings/illumina_midas31.sam

echo ""
echo "Mapped reads. Removing supplementary mappings based on SAM bit flags"
samtools view -F 256 -F 4 -F 2048 output/mappings/illumina_midas37.sam -o output/mappings/il_midas37_nodupes.sam
date

sed '/^@/ d' output/mappings/il_midas37_nodupes.sam | \
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
  }' > idmapped_ilv13.txt

awk -F"\t" -v OFS='\t' '{print $1, $3}' idmapped_ilv13.txt | sed 's/Zbarcode=/\t/' >  output/mappings/minimap_il_midas37_processed.txt
